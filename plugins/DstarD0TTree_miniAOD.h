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
		void initialize();



		// ----------member data ---------------------------
		bool doMC, doRec, debug, selectionCuts, triggerOn;
		double pi_mass, k_mass;
      std::string triggerName_;
		std::vector<int> dScandsKpi;
		std::vector<reco::TransientTrack>  goodTracks;
		std::vector<reco::TransientTrack>  goodTracksD0;
		std::vector<reco::TransientTrack> slowPiTracks;
		std::vector<reco::TransientTrack> t_tks;
		std::vector<reco::TransientTrack> tksD0;
		TTree *data;

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

		TransientVertex v_D0;
		math::XYZTLorentzVector d0kpi_p4;

		double Ebeam_;		double comEnergy_;
		double DstarSignificance3D_ ;
		double D0Significance3D_;

		//counters
		//int TotalTracks, TracksAfterTrigger, TracksHasTrackDetails, TracksChargeZero, TracksEta, TracksHighPurity, TracksPDG211;
		//int TracksPtZeroFive, TracksChi3, TracksNumberOfHits2, TracksDxyThree, TracksDzThree, TrackSlowPionCandidates;
		//int TracksPtZeroSix, TracksChiTwoFive, TracksNumberOfHits5, TracksNumberOfPixelHits2, TracksDxyZeroOne, TracksDzOne; 
		//int TrackKaonPionCandidates, D0AfterLorentzVector, DsMinusD0Zerothree, TransientTrackOfpiK, SVConfidenceLevel;
		//int D0AfterLorentzVectorKarman, PointingcosPhi, Significance, D0MinusPDG, D0pTThree, DsAfterLorentzVector, DsMinusD0, Observation;

		//counters
		unsigned long TotalTracks, TracksHasTrackDetails, TracksChargeZero, TracksEta, TracksHighPurity, TracksPDG211;
		unsigned long TracksPtZeroFive, TracksChi3, TracksNumberOfHits2, TracksDxyThree, TracksDzThree, TrackSlowPionCandidates;
		unsigned long TracksPtZeroSix, TracksChiTwoFive, TracksNumberOfHits5, TracksNumberOfPixelHits2, TracksDxyZeroOne, TracksDzOne; 
		unsigned long TrackKaonPionCandidates, D0AfterLorentzVector, D0MinusPDGzero2, DsMinusD0Zerothree, TransientTrackOfpiK, SVConfidenceLevel;
		unsigned long D0AfterLorentzVectorKarman, PointingcosPhi, Significance, D0MinusPDG, D0pTThree, DsAfterLorentzVector, DsMinusD0, Observation;
		unsigned long D0Candidates, DsCandidates;

		unsigned long TracksPt0p8, D0PromptTracksChi5, D0PromptTracksNumberHits5, D0PromptTracksPixelHits2;
		unsigned long D0PromptTracksDxy0p1, D0PromptTracksDz0p5, TracksD0combination;
		unsigned long D0PromptminusPDG1p0, D0PromptSVConfidenceLevel, D0PromptPointingCosPhi; 
		unsigned long D0PromptSignificance, D0PromptCandidates, D0PromptTracksEta2p5;
		unsigned long D0PromptKpiAfterTransientp0;

		int Total_Events,runNumber,eventNumber,lumi;
		int Total_Events_test, Triggered_Event, Triggered_Event_test;
		int ND0KpiCand, NKpiCand,NdsKpiMC,FlagMC,FlagRec,n_pVertex,ntracksD0Kpi,ntracksDstar,HLTPath_ ;   

		double PVx,PVy,PVz,PVerrx,PVerry,PVerrz,lumiWeight_;  
		double HFEnergyMinus,HFEnergyPlus; 
		TLorentzVector vD0_1, vD0_2, vD0kaon, vD0pion;
		double mass1, mass2;
		bool comb1, comb2, combOR;     
		// std::vector<int> FlagMC;

		std::vector<std::string> NameTrigger;
		//std::vector<TString> NameTrigger;
		std::vector<double> D0Kpi_VtxProb, D0Kpipt, D0Kpieta, D0Kpiphi, D0Kpi_VtxPosx, D0Kpi_VtxPosy, 
		D0Kpi_VtxPosz, D0Kpi_Vtxerrx, D0Kpi_Vtxerry, D0Kpi_Vtxerrz, D0Kpi_DispAngle;
		std::vector<double> D0Kpimass,TrkD0Keta,TrkD0pieta,TrkD0Kphi,TrkD0piphi;
		std::vector<double> TrkD0Kdxy, TrkD0pidxy, TrkD0Kdz, TrkD0pidz, TrkD0Knhits, TrkD0pinhits, 
		TrkD0Kchi2, TrkD0pichi2, D0DeltaR, TrkD0Kpt, TrkD0pipt, D0KpisXY_vec, D0Kpis3D_vec, D0_kT_vec;

		std::vector<double> D0_VtxProb, D0pt, Dspt, D0eta, Dseta, D0phi, Dsphi, D0_VtxPosx, D0_VtxPosy, D0_VtxPosz, 
		D0_Vtxerrx, D0_Vtxerry, D0_Vtxerrz, TrkKdxy;
		std::vector<double> Dsmass, Trkpidxy, TrkSdxy, TrkKdz, Trkpidz, TrkSdz, TrkKnhits, Trkpinhits, TrkSnhits, 
		TrkKchi2, Trkpichi2, TrkSchi2, DSDeltaR, TrkKpt,Trkpipt;
		std::vector<double> D0mass, TrkKmass, Trkpimass, TrkSmass, TrkSpt, TrkKeta, Trkpieta, TrkSeta, TrkKphi, 
		Trkpiphi, TrkSphi, TrkScharge, D0fromDSsXY_vec;
		std::vector<double> D0fromDSs3D_vec, Anglephi_vec;

		std::vector<double> D0_VtxProbWrong, D0ptWrong, DsptWrong, D0etaWrong, DsetaWrong, D0phiWrong, DsphiWrong, D0_VtxPosxWrong, D0_VtxPosyWrong, D0_VtxPoszWrong, 
		D0_VtxerrxWrong, D0_VtxerryWrong, D0_VtxerrzWrong, TrkKdxyWrong;
		std::vector<double> DsmassWrong, TrkpidxyWrong, TrkSdxyWrong, TrkKdzWrong, TrkpidzWrong, TrkSdzWrong, TrkKnhitsWrong, TrkpinhitsWrong, TrkSnhitsWrong, 
		TrkKchi2Wrong, Trkpichi2Wrong, TrkSchi2Wrong, DSDeltaRWrong, TrkKptWrong, TrkpiptWrong;
		std::vector<double> D0massWrong, TrkKmassWrong, TrkpimassWrong, TrkSmassWrong, TrkSptWrong, TrkKetaWrong, TrkpietaWrong, TrkSetaWrong, TrkKphiWrong, 
		TrkpiphiWrong, TrkSphiWrong, TrkSchargeWrong, D0fromDSsXY_vecWrong;
		std::vector<double> D0fromDSs3D_vecWrong, Anglephi_vecWrong;

		std::vector<double> MxFromPFCands_,EPlusPzFromPFCands_,EMinusPzFromPFCands_,sumEHFPlusFromPFCands_,sumEHFMinusFromPFCands_;
		std::vector<double> xiPlusFromPFCands_,xiMinusFromPFCands_,etaMaxFromPFCands_,etaMinFromPFCands_;
		std::vector<double> etaMaxFromPFCandsNew_,etaMinFromPFCandsNew_,missingMassFromXiFromPFCands_;

		std::vector<double> MCDseta,MCDsphi,MCDspt,MCDsenergy,MCDsp,MCDset,MCDsrapidity,MCDsmass;
		std::vector<double> MCD0eta,MCD0phi,MCD0pt,MCD0energy,MCD0p,MCD0et,MCD0rapidity,MCD0mass;

		std::vector<double> MCpromptD0eta, MCpromptD0phi, MCpromptD0pt, MCpromptD0energy, 
		MCpromptD0p, MCpromptD0et, MCpromptD0rapidity, MCpromptD0mass;
		std::vector<double> MCpromptD0_Keta, MCpromptD0_Kphi, MCpromptD0_Kpt, MCpromptD0_Kenergy,
		 MCpromptD0_Kp, MCpromptD0_Ket, MCpromptD0_Krapidity, MCpromptD0_Kmass;
		std::vector<double> MCpromptD0_Pieta, MCpromptD0_Piphi, MCpromptD0_Pipt, MCpromptD0_Pienergy, 
		MCpromptD0_Pip, MCpromptD0_Piet, MCpromptD0_Pirapidity, MCpromptD0_Pimass;
		std::vector<double> MCpromptD0_DispAngle,MCpromptD0_Kt;

		//      std::vector<double> xiGenPlus_,xiGenMinus_,MxGen_,MxGenRange_,sumEnergyHEPlusGen_,sumEnergyHEMinusGen_,sumEnergyHFPlusGen_;
		//    std::vector<double> sumEnergyHFMinusGen_,etaMaxGen_,etaMinGen_,deltaEtaGen_,etaGapLow_,etaGapHigh_,MxGenPlus_,MxGenMinus_;

		std::vector<double> pfsis1Eta_max,pfsis2Eta_max,pfsis1Eta_min,pfsis2Eta_min,deltaEtapf;

};
#endif 

