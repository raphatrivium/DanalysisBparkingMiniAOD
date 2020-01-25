#ifndef DstarD0TTree_DstarD0TTree_h
#define DstarD0TTree_DstarD0TTree_h


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

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
#include "Geometry/Records/interface/PTrackerParametersRcd.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector.h"

#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
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

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

//Weighting PileUp
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SMPJ/AnalysisFW/interface/QCDEventHdr.h"
#include "SMPJ/AnalysisFW/interface/QCDEvent.h"

//JetCollection (CaloTower) - Gapside
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

//-------------------------------------------------------------
#include "DataFormats/CaloTowers/interface/CaloTower.h"
//#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//--------------------------------------------------------------

//Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//#include "DstarD0/AnaDStarD0/interface/FWLiteTools.h"
//#include "ForwardAnalysis/Utilities/interface/LargestGenRapidityGap.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TLorentzVector.h"
//#include  "DStarD0/DStarD0Analysis/EventData.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <string>

//
// class declaration
//

class DstarD0TTree : public edm::EDAnalyzer {
	public:
		explicit DstarD0TTree(const edm::ParameterSet&);
		~DstarD0TTree();
		void RecDstar(const edm::Event& iEvent, const edm::EventSetup&, const reco::Vertex& RecVtx);
		void RecDstarWrongCombination(const edm::Event& iEvent, const edm::EventSetup&, const reco::Vertex& RecVtx);
		void RecD0(const edm::Event& iEvent, const edm::EventSetup&, const reco::Vertex& RecVtx);
		double FindAngle(const reco::Vertex& , const TransientVertex& , const math::XYZTLorentzVector& ) ;
		//double FindAngleMCpromptD0(const GenParticle&);


	private:

		virtual void beginJob() ;
		virtual void endJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		void GenDstarInfo(const edm::Event& iEvent,const edm::EventSetup&);
		void GenD0Info(const edm::Event& iEvent,const edm::EventSetup&);

		void assignStableDaughters(const reco::Candidate* p, std::vector<int> & pids);
		bool TriggerInfo(const edm::Event&, edm::Handle<edm::TriggerResults>, edm::Handle<pat::PackedTriggerPrescales> ,TString trigname);
		std::vector<std::string> TriggersFired(const edm::Event&, edm::Handle<edm::TriggerResults>, edm::Handle<pat::PackedTriggerPrescales>);
		void initialize();

		// ----------member data ---------------------------
		bool doMC, doRec, debug, selectionCuts, triggerOn, triggerReferenceOn, TracksOn;
		double pi_mass, k_mass;
      std::string triggerName_;
		std::string triggerReferenceName_;
		std::string triggerReferenceName2_;
		std::vector<int> dScandsKpi;
		std::vector<reco::TransientTrack>  goodTracks;
		std::vector<reco::TransientTrack>  goodTracksD0;
		std::vector<reco::TransientTrack> slowPiTracks;
		//std::vector<reco::TransientTrack> t_tks;
		//std::vector<reco::TransientTrack> tksD0;

		TTree *data;

		//Weighting
		edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSumInfoInputTag_;
		edm::LumiReWeighting *LumiWeights_;
		edm::EDGetTokenT<GenEventInfoProduct> mEventInfo_;

		

		//Triggers
		edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
		edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
		edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescales_;

		//MiniAOD
		edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trkToken_; 
		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;

		edm::EDGetTokenT<reco::GenParticleCollection> genParticlesTokenDstar_;
		edm::EDGetTokenT<reco::GenParticleCollection> genParticlesTokenD0_;

		
		//edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
		//edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;

		//std::vector<TString> NamesFiredTrigger;
		//std::vector<std::vector<TString>> NameOfFiredTriggers;
		std::vector<std::string> NameOfFiredTriggers;

		QCDEvent *mEvent;

		double PUWeight;
		double nPU;
	
		double Ebeam_, comEnergy_, DstarSignificance3D_ , D0Significance3D_;

		int 	Total_Events, Triggered_Event, runNumber,eventNumber,lumi, Total_Events_test,
				Triggered_Event_test, TriggeredReference_Event_test, ND0KpiCand, NKpiCand, 
				NdsKpiMC, FlagMC, FlagRec, n_pVertex, ntracksD0Kpi, ntracksDstar,
				TriggeredReference_Event_test2;  

		double PVx,PVy,PVz,PVerrx,PVerry,PVerrz,lumiWeight_,  
				 HFEnergyMinus, HFEnergyPlus, mass1, mass2;

		//counters
		//int TotalTracks, TracksAfterTrigger, TracksHasTrackDetails, TracksChargeZero, TracksEta, TracksHighPurity, TracksPDG211;
		//int TracksPtZeroFive, TracksChi3, TracksNumberOfHits2, TracksDxyThree, TracksDzThree, TrackSlowPionCandidates;
		//int TracksPtZeroSix, TracksChiTwoFive, TracksNumberOfHits5, TracksNumberOfPixelHits2, TracksDxyZeroOne, TracksDzOne; 
		//int TrackKaonPionCandidates, D0AfterLorentzVector, DsMinusD0Zerothree, TransientTrackOfpiK, SVConfidenceLevel;
		//int D0AfterLorentzVectorKarman, PointingcosPhi, Significance, D0MinusPDG, D0pTThree, DsAfterLorentzVector, DsMinusD0, Observation;

		//----------------------------
		//Tracks
		//----------------------------
		unsigned long 	CounterTotalTracks, CounterTracksHasTrackDetails, CounterTracksCharged, 
							CounterTracksEta, CounterTracksHighPurity, CounterTracksPDG211,
		 					CounterTracksPtZeroFive, CounterTracksChi3, CounterTracksNumberOfHits2, 
							CounterTracksDxyThree, CounterTracksDzThree, CounterTrackSlowPionCandidates,
		 					CounterTracksPtZeroSix, CounterTracksChiTwoFive, CounterTracksNumberOfHits5, 
							CounterTracksNumberOfPixelHits2, CounterTracksDxyZeroOne, CounterTracksDzOne,
							CounterTrackKaonPionCandidates;

		std::vector<double> 	TracksCharge, TracksEta, TracksPt, TracksChi2, TracksNumberOfHits, TracksPhi,
									TracksdxyError, TracksdzError, TracksNumberOfPixelHits, Tracksdxy, Tracksdz;

		//----------------------------
		//D* QUANTITIES
		//----------------------------
		unsigned long	CounterD0AfterLorentzVector, CounterD0MinusPDGzero2, CounterDsMinusD0Zerothree, 
							CounterTransientTrackOfpiK, CounterSVConfidenceLevel,	CounterD0AfterLorentzVectorKarman, 
							CounterPointingcosPhi, CounterSignificance3D, CounterD0MinusPDG, CounterD0pTThree, 
							CounterDsAfterLorentzVector, CounterDsMinusD0, CounterD0Candidates, 
							CounterDsCandidates;
		
		std::vector<double> 	D0Kpi_VtxProb, D0Kpipt, D0Kpieta, D0Kpiphi, D0Kpi_VtxPosx, D0Kpi_VtxPosy, 
									D0Kpi_VtxPosz, D0Kpi_Vtxerrx, D0Kpi_Vtxerry, D0Kpi_Vtxerrz, D0Kpi_DispAngle,
									D0Kpimass, TrkD0Keta, TrkD0pieta, TrkD0Kphi, TrkD0piphi, TrkD0Kdxy, TrkD0pidxy, 
									TrkD0Kdz, TrkD0pidz, TrkD0Knhits, TrkD0pinhits,	TrkD0Kchi2, TrkD0pichi2, 
									D0DeltaR, TrkD0Kpt, TrkD0pipt, D0KpisXY_vec, D0Kpis3D_vec, D0_kT_vec;

		std::vector<double> 	D0_VtxProb, D0pt, Dspt, D0eta, Dseta, D0phi, Dsphi, D0_VtxPosx, D0_VtxPosy,
									D0_VtxPosz, D0_Vtxerrx, D0_Vtxerry, D0_Vtxerrz, TrkKdxy, Dsmass, Trkpidxy, 
									TrkSdxy, TrkKdz, Trkpidz, TrkSdz, TrkKnhits, Trkpinhits, TrkSnhits, 
									TrkKchi2, Trkpichi2, TrkSchi2, DSDeltaR, TrkKpt,Trkpipt,	D0mass, TrkKmass, 
									Trkpimass, TrkSmass, TrkSpt, TrkKeta, Trkpieta, TrkSeta, TrkKphi, 
									Trkpiphi, TrkSphi, TrkScharge, D0fromDSsXY_vec,	D0fromDSs3D_vec, Anglephi_vec, 
									D0fromDSd3D_vec, D0fromDSe3D_vec, D0Kpid3D_vec, D0Kpie3D_vec,
									D0fromDSdXY_vec, D0fromDSeXY_vec, D0KpidXY_vec, D0KpieXY_vec; 
									

		//----------------------------
		//D* QUANTITIES WRONG COMBINATION
		//----------------------------
		std::vector<double> 	D0_VtxProbWrong, D0ptWrong, DsptWrong, D0etaWrong, DsetaWrong, D0phiWrong, 
									DsphiWrong, D0_VtxPosxWrong, D0_VtxPosyWrong, D0_VtxPoszWrong,	D0_VtxerrxWrong, 
									D0_VtxerryWrong, D0_VtxerrzWrong, TrkKdxyWrong,	DsmassWrong, TrkpidxyWrong, 
									TrkSdxyWrong, TrkKdzWrong, TrkpidzWrong, TrkSdzWrong, TrkKnhitsWrong, 
									TrkpinhitsWrong, TrkSnhitsWrong, TrkKchi2Wrong, Trkpichi2Wrong, TrkSchi2Wrong, 
									DSDeltaRWrong, TrkKptWrong, TrkpiptWrong, D0massWrong, TrkKmassWrong, 
									TrkpimassWrong, TrkSmassWrong, TrkSptWrong, TrkKetaWrong, TrkpietaWrong, 
									TrkSetaWrong, TrkKphiWrong, TrkpiphiWrong, TrkSphiWrong, TrkSchargeWrong, 
									D0fromDSsXY_vecWrong, D0fromDSs3D_vecWrong, Anglephi_vecWrong;

		//----------------------------
		//D* MC
		//----------------------------
		std::vector<double> 	MCDseta,MCDsphi,MCDspt,MCDsenergy,MCDsp,MCDset,MCDsrapidity,MCDsmass,
		 							MCD0eta,MCD0phi,MCD0pt,MCD0energy,MCD0p,MCD0et,MCD0rapidity,MCD0mass,
							 		MCDsKeta,MCDsKphi,MCDsKpt,MCDsKenergy,MCDsKp,MCDsKet,MCDsKrapidity,
									MCDsKmass, MCDsPieta,MCDsPiphi,MCDsPipt, MCDsPienergy, MCDsPip, 
									MCDsPiet, MCDsPirapidity, MCDsPimass, Dseta_vec, MCDseta_vec, 
									Dsphi_vec, MCDsphi_vec, Dspt_vec, MCDspt_vec, D0fromDsmass_vec,
									deltaRDs_vec;

		//----------------------------
		//D0 QUANTITIES
		//----------------------------
		unsigned long 	CounterTracksPt0p8, CounterD0PromptTracksChi5, CounterD0PromptTracksNumberHits5, 
							CounterD0PromptTracksPixelHits2,CounterD0PromptTracksDxy0p1, CounterD0PromptTracksDz0p5,
							CounterTracksD0combination, CounterD0PromptminusPDG1p0, CounterD0PromptSVConfidenceLevel,
							CounterD0PromptPointingcosPhi, CounterD0PromptSignificance3D, CounterD0PromptCandidates, 
							CounterD0PromptTracksEta2p5, CounterD0PromptKpiAfterTransientp0;

		bool comb1, comb2;
		

		//----------------------------
		//D0 MC
		//----------------------------			
		std::vector<double> 	MCpromptD0eta, MCpromptD0phi, MCpromptD0pt, MCpromptD0energy, 
									MCpromptD0p, MCpromptD0et, MCpromptD0rapidity, MCpromptD0mass,
		 							MCpromptD0_Keta, MCpromptD0_Kphi, MCpromptD0_Kpt, MCpromptD0_Kenergy,
		 							MCpromptD0_Kp, MCpromptD0_Ket, MCpromptD0_Krapidity, MCpromptD0_Kmass,
		 							MCpromptD0_Pieta, MCpromptD0_Piphi, MCpromptD0_Pipt, MCpromptD0_Pienergy, 
									MCpromptD0_Pip, MCpromptD0_Piet, MCpromptD0_Pirapidity, MCpromptD0_Pimass,
		 							MCpromptD0_DispAngle,MCpromptD0_Kt;

		std::vector<double> 	D0eta_vec, MCD0eta_vec, D0phi_vec, MCD0phi_vec, D0pt_vec, 
									MCD0pt_vec, D0Kt_vec, D0Sxy_vec, D0OpAngle_vec, deltaRD0_vec;
	
		

	
};
#endif 

