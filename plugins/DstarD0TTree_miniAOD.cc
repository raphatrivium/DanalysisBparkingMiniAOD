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
	triggerName_(iConfig.getUntrackedParameter<std::string>("PathName","HLT_Mu9_IP6_part0_v2")),        
 	//Triggers
	triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))),
	triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
	trkToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("tracks"))), //MiniAOD
	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("recVtxs"))),
	genParticlesTokenDstar_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gens"))),
	genParticlesTokenD0_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gensD0"))),
	comEnergy_(iConfig.getParameter<double>("comEnergy"))
{     
	counter = 0;
	countInTriggered = 0;
	Ebeam_ = comEnergy_/2.;
	edm::Service<TFileService> fs;
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
    counter++; 
	//To clear and initialize variables
	initialize();

	//run, event, lumi section
	runNumber= iEvent.id().run();
	eventNumber= iEvent.id().event();
	lumi= iEvent.luminosityBlock();
	
	Handle<double> lumiWeight;
	iEvent.getByLabel("lumiWeight",lumiWeight);

	lumiWeight_=lumi;
	
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
   if(triggerFired) countInTriggered++;

	NameTrigger.push_back(triggerName_);
	
	// Getting tracks from vert(ex)ices
	edm::Handle< View < pat::PackedCandidate >> tracks; //access that PackedCandidate
	iEvent.getByToken(trkToken_,tracks);

	//miniAOD
	edm::ESHandle<TransientTrackBuilder> theB; 
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

    //Only events in which the path actually fired had stored the filter results and products:	  
    //bool triggerFired = TriggerInfo(iEvent,triggerBits,triggerPrescales,triggerName_);
    //if(triggerFired) countInTriggered++;

	//Loop in the trakcs of the PackedCandidate 
	for(View<pat::PackedCandidate>::const_iterator iTrack1 = tracks->begin(); iTrack1 != tracks->end(); ++iTrack1 ) 
	{   TotalTracks++;
        if (triggerFired) 
			{ 	TracksAfterTrigger++;
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
            }
			
				Observation++;
				//Kaon and Pion Candidates
				if( iTrack1->pt()>0.5 /*&& (fabs(iTrack1->pdgId()) == 211)*/)
				{
                  TracksPtZeroSix++;
                  if(iTrack1->pseudoTrack().normalizedChi2() > 2.5) continue;
                  //cout << "iTrack1->pseudoTrack().normalizedChi2(): "<< iTrack1->pseudoTrack().normalizedChi2() << endl ;
                  //cout << "iTrack1->pseudoTrack().Chi2(): "<< iTrack1->pseudoTrack().chi2() << endl ;
                  TracksChiTwoFive++;
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
				}

            // SELECTING TRACKS FOR D0             
				if(fabs(iTrack1->eta())<2.5 && iTrack1->pseudoTrack().normalizedChi2() < 5.0 && iTrack1->pseudoTrack().hitPattern().numberOfValidHits() >= 5 && iTrack1->pseudoTrack().hitPattern().numberOfValidPixelHits() >= 2 && iTrack1->pt() > 0.5 && iTrack1->p() >1.0 && fabs(iTrack1->dz())<0.5 && fabs(iTrack1->dxy())<0.1)
        		{
					reco::TransientTrack  D0TT = theB->build(iTrack1->pseudoTrack());
					goodTracksD0.push_back(D0TT);
            }
			}//trigger
   }//loop packed candidates
//}//trigger

	if (debug) cout << " goodTracks size " << goodTracks.size() << endl;
    ntracksDstar = slowPiTracks.size();
    ntracksD0Kpi = goodTracksD0.size();

	//Vertex Informations
   if(triggerFired)
   {
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
      RecD0(iEvent,iSetup,RecVtx);    //Reconstruction of prompt D0
      GenDstarInfo(iEvent,iSetup); //Stores information from D0 and its products.
      GenD0Info(iEvent,iSetup); //Stores information from D0 and its products.

		//FindAngle(RecVtx,v_D0,d0kpi_p4); //Calculates the opening angle
      //FindAngleMCpromptD0(p);
		data->Fill();
    }
}

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
				
				TransientTrack K;
				TransientTrack pi;
                                
				//Right Combination Charges
				if(trk1.charge() == trkS.charge())
				{ //  cout << "Right Combination Charges" << endl;
					pi = trk1;
					K = trk2;
				}
				else
				{
					K = trk1;
					pi = trk2;
				}

				//Wrong Combination Charges - background
				/*if(trk1.charge() == trkS.charge())
				{ //  cout << "Right Combination Charges" << endl;
					pi = trk2;
					K = trk1;
				}
				else
				{
					K = trk2;
					pi = trk1;
				}	*/	
                               
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

				//cout << "-------------------------------------------" << endl;
				//cout << "2 Trk S chi2: " << trkS.track().normalizedChi2() << endl;
				//cout << "2 Trk pi chi2: " << pi.track().normalizedChi2() << endl;
				//cout << "2 Trk K chi2: " << K.track().normalizedChi2() << endl;
				//cout << "2 Trk pi_f chi2: " << pi_f.track().normalizedChi2() << endl;
				//cout << "2 Trk K_f normalizedChi2: " << K_f.track().normalizedChi2() << endl;
				//cout << "2 Trk K_f chi2: " << K_f.track().chi2() << endl;
				
				//SV Confidence Level
				if(vtxProb < 0.01) continue;
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
				
  				double anglephi = FindAngle(RecVtx,v,d0_p4);
				double cosPhi = cos(anglephi);
				if( cosPhi < 0.99 ) continue;
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

				TrkKpt.push_back(K_f.track().pt());
				Trkpipt.push_back(pi_f.track().pt());
				TrkSpt.push_back(trkS.track().pt());

				//cout << "significance 3D: "<< D0fromDSs3D << endl;
				//if( D0fromDSs3D < 3. ) continue;
				if( D0fromDSs3D < .3 ) continue;
				Significance++;

				//Difference between D0 and D0PDG 5sigmas
				if( fabs(d0mass - 1.86484) > 0.1 ) continue;
				D0MinusPDG++;

				//Pt Cut of mesnons D0
				if( d0_p4.Pt() < 3. ) continue;
				D0pTThree++;
								
				D0Candidates++;

				math::XYZTLorentzVector dS_p4 = d0_p4 + p4_S;
				double dsmass = dS_p4.M();
				DsAfterLorentzVector++;
								                      											
				//Difference between D* and D0 -> Must be close to pion
				if( (dsmass - d0mass) > 0.16) continue;
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
									
				if(NKpiCand>999) break;
			}
			if(NKpiCand>999) break;
	 } 
	if(NKpiCand>999) break;
	}
}//End RecDstar

//***********************************************************************************
void DstarD0TTree::assignStableDaughters(const reco::Candidate* p, std::vector<int> & pids){
	//Evaluate if a particle is stable (status 1). If not it repeat the process until find the a stable one.
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
		TransientTrack trk1D0 = goodTracks[i];

		for(size_t j=i+1;j<goodTracksD0.size();j++)
		{
			TransientTrack trk2D0 = goodTracksD0[i];
			//Testing charge and if tracks are equal
			if(trk1D0 == trk2D0) continue;
			if(trk1D0.charge() == trk2D0.charge()) continue;

			//Recosntruction of D0 momentum
			math::XYZVector D0_p = trk1D0.track().momentum() + trk2D0.track().momentum();

			//TransientTrack KfromD0 = 0, PifromD0 = 0;	
			TransientTrack KfromD0 , PifromD0 ;	        
            
			comb1 = comb2 = false ;
			mass1 = mass2 = 0 ;
			combOR = false;

			//Reconstructio of 4-vector D0 Prompt
			vD0kaon.SetPtEtaPhiM(trk1D0.track().pt(), trk1D0.track().eta(), trk1D0.track().phi(), k_mass);
			vD0pion.SetPtEtaPhiM(trk2D0.track().pt(), trk2D0.track().eta(), trk2D0.track().phi(), pi_mass);
			vD0_1 = vD0kaon + vD0kaon;
			mass1 = vD0_1.M();

			//Reconstruction of 4-vector D0 Prompt
			vD0pion.SetPtEtaPhiM(trk1D0.track().pt(), trk1D0.track().eta(), trk1D0.track().phi(), pi_mass);
			vD0kaon.SetPtEtaPhiM(trk2D0.track().pt(), trk2D0.track().eta(), trk2D0.track().phi(), k_mass);
			vD0_2 = vD0kaon + vD0kaon;
			mass2 = vD0_2.M();

			//cout << "mass1: " << mass1 << "------- mass2: " <<  mass2 << endl;

			if( fabs(mass1-1.86484) < 1.0) comb1 = true;		
			if( fabs(mass2-1.86484) < 1.0)  comb2 = true;

			//====================================
			//If one of the conditions is satisfied	
			if(comb1 || comb2)
			{
				if (comb1)
				{
					//cout << "comb1: " << endl;
					KfromD0 = trk1D0; PifromD0 = trk2D0;
					math::XYZVector K_p = trk1D0.track().momentum();

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

					math::XYZTLorentzVector p4_KfromD0(KfromD0_f.track().px(),
																	KfromD0_f.track().py(),
																	KfromD0_f.track().pz(),
																	sqrt(pow(KfromD0_f.track().p(),2)+pow(k_mass,2))  );
					math::XYZTLorentzVector p4_pifromD0(pifromD0_f.track().px(),
																	pifromD0_f.track().py(),
																	pifromD0_f.track().pz(),
																	sqrt(pow(pifromD0_f.track().p(),2)+pow(pi_mass,2)) );  
					d0kpi_p4 = p4_KfromD0 + p4_pifromD0;

					double D0_kT = sqrt(  (K_p).Cross(D0_p).Mag2() / D0_p .Mag2() ) ;
					double d0kpimass = d0kpi_p4.M();
					//cout << "d0kpimass comb1: " << d0kpimass << endl;
					if(fabs(d0kpimass - 1.86484)>0.15) continue;
					ND0KpiCand++;

					double dispAngle = FindAngle(RecVtx,v_D0,d0kpi_p4);

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

				} // end of combination1

				// -----
				// end of comb1

				if(comb2)
				{
					//cout << "comb2: " << endl;
					KfromD0 = trk2D0; PifromD0 = trk1D0;
					math::XYZVector K_p = trk2D0.track().momentum();

					vector<TransientTrack> tksD0;
					tksD0.push_back(KfromD0);
					tksD0.push_back(PifromD0);
					KalmanVertexFitter kalman(true);
					v_D0 = kalman.vertex(tksD0);

					if(!v_D0.isValid() || !v_D0.hasRefittedTracks()) continue;
					double D0KpivtxProb =TMath::Prob( (Double_t) v_D0.totalChiSquared(), (Int_t) v_D0.degreesOfFreedom());
					TransientTrack KfromD0_f = v_D0.refittedTrack(KfromD0);
					TransientTrack pifromD0_f = v_D0.refittedTrack(PifromD0);      
					VertexDistanceXY vD0KpidXY ;			
					double D0KpidXY = vD0KpidXY.distance(RecVtx,v_D0).value() ;
					double D0KpieXY = vD0KpidXY.distance(RecVtx,v_D0).error() ;
					double D0KpisXY = D0KpidXY / D0KpieXY;
					VertexDistance3D vD0Kpid3D ;
					double D0Kpid3D = vD0Kpid3D.distance(RecVtx,v_D0).value() ;
					double D0Kpie3D = vD0Kpid3D.distance(RecVtx,v_D0).error() ;
					double D0Kpis3D = D0Kpid3D / D0Kpie3D;

					math::XYZTLorentzVector p4_KfromD0(KfromD0_f.track().px(),
																	KfromD0_f.track().py(),
																	KfromD0_f.track().pz(),
																	sqrt(pow(KfromD0_f.track().p(),2)+pow(k_mass,2)));
					math::XYZTLorentzVector p4_pifromD0(pifromD0_f.track().px(),
																	pifromD0_f.track().py(),
																	pifromD0_f.track().pz(),
																	sqrt(pow(pifromD0_f.track().p(),2)+pow(pi_mass,2)));  
					d0kpi_p4 = p4_KfromD0 + p4_pifromD0;
					double D0_kT = sqrt(  (K_p).Cross(D0_p).Mag2() / D0_p .Mag2() ) ;
					double d0kpimass = d0kpi_p4.M();
					//cout << "d0kpimass comb2: " << d0kpimass << endl;
					if(fabs(d0kpimass - 1.86484)>0.15) continue;
					ND0KpiCand++;

					double dispAngle = FindAngle(RecVtx,v_D0,d0kpi_p4);

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

				}// end of combination 2

				combOR = true;
			} //end of comb1 OR comb 2

			//=============================
			//If both meet the requiremets			
			if(comb1 && comb2 &! combOR){
				//cout << "comb 1 and 2: " << endl;

				if(fabs( mass1-1.864 ) < fabs( mass2-1.864 ))
				{
					KfromD0 = trk1D0; PifromD0 = trk2D0;

					math::XYZVector K_p = trk1D0.track().momentum();

					vector<TransientTrack> tksD0;
					tksD0.push_back(KfromD0);
					tksD0.push_back(PifromD0);
					KalmanVertexFitter kalman(true);
					v_D0 = kalman.vertex(tksD0);

					if(!v_D0.isValid() || !v_D0.hasRefittedTracks()) continue;
					double D0KpivtxProb =TMath::Prob( (Double_t) v_D0.totalChiSquared(), (Int_t) v_D0.degreesOfFreedom());
					TransientTrack KfromD0_f = v_D0.refittedTrack(KfromD0);
					TransientTrack pifromD0_f = v_D0.refittedTrack(PifromD0);
					VertexDistanceXY vD0KpidXY ;
					double D0KpidXY = vD0KpidXY.distance(RecVtx,v_D0).value() ;
					double D0KpieXY = vD0KpidXY.distance(RecVtx,v_D0).error() ;
					double D0KpisXY = D0KpidXY / D0KpieXY;
					VertexDistance3D vD0Kpid3D ;
					double D0Kpid3D = vD0Kpid3D.distance(RecVtx,v_D0).value() ;
					double D0Kpie3D = vD0Kpid3D.distance(RecVtx,v_D0).error() ;
					double D0Kpis3D = D0Kpid3D / D0Kpie3D;

					math::XYZTLorentzVector p4_KfromD0(KfromD0_f.track().px(),
																	KfromD0_f.track().py(),
																	KfromD0_f.track().pz(),
																	sqrt(pow(KfromD0_f.track().p(),2)+pow(k_mass,2)));
					math::XYZTLorentzVector p4_pifromD0(pifromD0_f.track().px(),
																	pifromD0_f.track().py(),
																	pifromD0_f.track().pz(),
																	sqrt(pow(pifromD0_f.track().p(),2)+pow(pi_mass,2)));
					d0kpi_p4 = p4_KfromD0 + p4_pifromD0;
					double D0_kT = sqrt(  (K_p).Cross(D0_p).Mag2() / D0_p .Mag2() ) ;
					double d0kpimass = d0kpi_p4.M();
					//cout << "d0kpimass comb1_2: " << d0kpimass << endl;
					if(fabs(d0kpimass - 1.86484)>0.15) continue;
					ND0KpiCand++;

					double dispAngle = FindAngle(RecVtx,v_D0,d0kpi_p4);
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

				}

				if(fabs( mass2-1.864 ) < fabs( mass1-1.864 ))
				{
					KfromD0 = trk2D0; PifromD0 = trk1D0;

					math::XYZVector K_p = trk2D0.track().momentum();

					vector<TransientTrack> tksD0;
					tksD0.push_back(KfromD0);
					tksD0.push_back(PifromD0);
					KalmanVertexFitter kalman(true);
					v_D0 = kalman.vertex(tksD0);

					if(!v_D0.isValid() || !v_D0.hasRefittedTracks()) continue;
					double D0KpivtxProb =TMath::Prob( (Double_t) v_D0.totalChiSquared(), (Int_t) v_D0.degreesOfFreedom());
					TransientTrack KfromD0_f = v_D0.refittedTrack(KfromD0);
					TransientTrack pifromD0_f = v_D0.refittedTrack(PifromD0);
					VertexDistanceXY vD0KpidXY ;
					double D0KpidXY = vD0KpidXY.distance(RecVtx,v_D0).value() ;
					double D0KpieXY = vD0KpidXY.distance(RecVtx,v_D0).error() ;
					double D0KpisXY = D0KpidXY / D0KpieXY;
					VertexDistance3D vD0Kpid3D ;
					double D0Kpid3D = vD0Kpid3D.distance(RecVtx,v_D0).value() ;
					double D0Kpie3D = vD0Kpid3D.distance(RecVtx,v_D0).error() ;
					double D0Kpis3D = D0Kpid3D / D0Kpie3D;

					math::XYZTLorentzVector p4_KfromD0(	KfromD0_f.track().px(),
																	KfromD0_f.track().py(),
																	KfromD0_f.track().pz(),
																	sqrt(pow(KfromD0_f.track().p(),2)+pow(k_mass,2)));
					math::XYZTLorentzVector p4_pifromD0(pifromD0_f.track().px(),
																	pifromD0_f.track().py(),
																	pifromD0_f.track().pz(),
																	sqrt(pow(pifromD0_f.track().p(),2)+pow(pi_mass,2)));
					d0kpi_p4 = p4_KfromD0 + p4_pifromD0;
					double D0_kT = sqrt( (K_p).Cross(D0_p).Mag2() / D0_p .Mag2() ) ;
					double d0kpimass = d0kpi_p4.M();
					//cout << "d0kpimass comb1: " << d0kpimass << endl;
					if(fabs(d0kpimass - 1.86484)>0.15) continue;
					ND0KpiCand++;

					double dispAngle = FindAngle(RecVtx,v_D0,d0kpi_p4);
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

				}
			} //end of comb1 AND comb 2		
		}
	}
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
	Total_Events = 0; TotalTracks = 0; TracksAfterTrigger = 0; TracksHasTrackDetails =0; TracksChargeZero = 0; TracksEta = 0; 
	TracksHighPurity = 0; TracksPDG211 = 0; TracksPtZeroFive = 0; TracksChi3 = 0; TracksNumberOfHits2 = 0; 
	TracksDxyThree = 0; TracksDzThree = 0; TrackSlowPionCandidates = 0; TracksPtZeroSix = 0; TracksChiTwoFive = 0;
	TracksNumberOfHits5 = 0; TracksNumberOfPixelHits2 = 0; TracksDxyZeroOne = 0; TracksDzOne = 0; TrackKaonPionCandidates = 0;
	Observation = 0;

	D0AfterLorentzVector = 0; D0MinusPDGzero2 = 0; DsMinusD0Zerothree = 0; TransientTrackOfpiK = 0; SVConfidenceLevel = 0; D0AfterLorentzVectorKarman = 0; 
	PointingcosPhi = 0;	Significance = 0; D0pTThree = 0; DsAfterLorentzVector = 0; D0MinusPDG = 0; DsMinusD0 = 0;

	signalpf = false;
	TTBit_32 = 0;
	TTBit_33 = 0;
	TTBit_34 = 0;

	TTBit_8 = 0;
	TTBit_9 = 0;
	TTBit_10 = 0;

	nHFPlus = 0; nHFMinus = 0;

	
	NameTrigger.clear();
	D0_VtxProb.clear(); D0mass.clear(); TrkKmass.clear(); Trkpimass.clear(); TrkSmass.clear();Dsmass.clear();D0pt.clear();Dspt.clear();D0eta.clear();
	D0phi.clear();Dseta.clear();Dsphi.clear();D0_VtxPosx.clear();D0_VtxPosy.clear();D0_VtxPosz.clear();D0_Vtxerrx.clear(); D0_Vtxerry.clear();
	D0_Vtxerrz.clear();TrkKdxy.clear();Trkpidxy.clear();
	TrkSdxy.clear();TrkKdz.clear();Trkpidz.clear();TrkSdz.clear();
	TrkKnhits.clear();Trkpinhits.clear();TrkSnhits.clear();
	TrkKchi2.clear();Trkpichi2.clear();TrkSchi2.clear();
	DSDeltaR.clear();TrkKpt.clear();Trkpipt.clear();TrkSpt.clear();
	TrkKeta.clear();Trkpieta.clear();TrkSeta.clear();TrkKphi.clear();Trkpiphi.clear();TrkSphi.clear();TrkScharge.clear();
	D0fromDSsXY_vec.clear(); D0fromDSs3D_vec.clear(); Anglephi_vec.clear();

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
    cout << "Total # events: " << counter << "Total # events triggered by " << triggerName_ << ": " << countInTriggered << endl;
    
}

//***********************************************************************************
void DstarD0TTree::beginJob(){
        
    //cout << "DsCandidates: " << DsCandidates << " NKpiCand: " << NKpiCand << endl;
 
	data->Branch("Total_Events",&Total_Events,"Total_Events/I"); 
	data->Branch("countInTriggered",&countInTriggered,"countInTriggered/I");
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

	data->Branch("procId",&procId,"procId/I");

	data->Branch("TotalTracks",&TotalTracks,"TotalTracks/L");
	data->Branch("TracksAfterTrigger",&TracksAfterTrigger,"TracksAfterTrigger/L");
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
	// D0 Variables
	//======================================================

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
