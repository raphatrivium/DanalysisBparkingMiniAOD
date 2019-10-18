//Program to extract data (variables and vector format) of the root file and create histograms. This root file is made of others root files, therefore, has multiple entries which we access with a diferent method than a normal root file.

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atof */
#include <math.h>       /* sin */
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <vector>
#ifndef ROOT_TLatex
#define ROOT_TLatex
#ifndef ROOTiosfwd
#include "Riosfwd.h"
#endif
#ifndef ROOT_TText
#include "TText.h"
#endif
#ifndef ROOT_TAttLine
#include "TAttLine.h"
#endif
//#include "RooCBShape.h"  //Crystal Ball
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
//#include "RooRealVar.h"
//#include "RooDataSet.h"
//#include "RooGaussian.h"
#include "TCanvas.h"
//#include "RooPlot.h"
#include "TAxis.h"
//using namespace RooFit ;


//gSystem->AddIncludePath("/home/raphael/Desktop/analysis_2018");

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int analysisB2019()
{
	//call a file for a histogram style (Optional)
	gROOT->LoadMacro("styleTDR.C"); 
	setTDRStyle();
	
	//Lorentz Vector
	//TLorentzVector mu_1;
	//TLorentzVector mu_2;

	//To get the counters (variables int)
	bool Conters = true ;
	//To get the vector and make histograms
	bool VectorsHistograms = false ;
	//It is bind to bool VectorsHistograms. False = Make only Mesons D histograms, True = Make allhistogrmas
	bool Meson = false;

			
	//Variables	and Vectors
	int Total_Events = 0;

	int TotalTracks = 0;
	int TracksAfterTrigger = 0;
	int TracksHasTrackDetails = 0;
	int TracksChargeZero = 0;
	int TracksHighPurity = 0;
	int TracksPDG211 = 0;
	int TracksPtZeroFive = 0;
	int TracksNumberOfHits2 = 0;
	int TracksDxyThree = 0;
	int TracksDzThree = 0;
	int TrackSlowPionCandidates = 0;
	int	Observation = 0;
	int TracksPtZeroSix = 0;
	int TracksChiTwoFive = 0;
	int TracksNumberOfHits5 = 0;
	int TracksNumberOfPixelHits2 = 0;
	int TracksDzOne = 0;
	int TracksDxyZeroOne = 0;
	int TrackKaonPionCandidates = 0;

	int D0AfterLorentzVector = 0;
	int DsMinusD0Zerothree = 0;
	int TransientTrackOfpiK = 0;
	int PointingcosPhi = 0;
	int Significance = 0;
	int D0pTThree = 0;
	int D0Candidates = 0;
	int DsAfterLorentzVector = 0;
	int DsMinusD0 = 0;
	int DsCandidates = 0;
	
	std::vector<double>* D0mass = 0.;
	std::vector<double>* Dsmass = 0.;
	std::vector<double>* D0_VtxProb = 0.;
	std::vector<double>* D0_VtxPosx = 0.;
	std::vector<double>* D0_VtxPosy = 0.;
	std::vector<double>* D0_VtxPosz = 0.;
	std::vector<double>* D0_Vtxerrx = 0.;
	std::vector<double>* D0_Vtxerry = 0.;
	std::vector<double>* D0_Vtxerrz = 0.;
	std::vector<double>* D0eta = 0.;
	std::vector<double>* D0phi = 0.;
	std::vector<double>* Dseta = 0.;
	std::vector<double>* Dsphi = 0.;
	std::vector<double>* TrkKpt = 0.;
	std::vector<double>* Trkpipt = 0.;
	std::vector<double>* TrkSpt = 0.;
	std::vector<double>* D0pt = 0.;
	std::vector<double>* Dspt = 0.;
	std::vector<double>* DSDeltaR = 0.;
	std::vector<double>* TrkKnhits = 0.;
	std::vector<double>* Trkpinhits = 0.;
	std::vector<double>* TrkSnhits = 0.;
	std::vector<double>* TrkKchi2 = 0.;
	std::vector<double>* Trkpichi2 = 0.;
	std::vector<double>* TrkSchi2 = 0.;
	std::vector<double>* TrkKdxy = 0.;
	std::vector<double>* Trkpidxy = 0.;
	std::vector<double>* TrkSdxy = 0.;
	std::vector<double>* TrkKdz = 0.;
	std::vector<double>* Trkpidz = 0.;
	std::vector<double>* TrkSdz = 0.;
	std::vector<double>* TrkKeta = 0.;
	std::vector<double>* Trkpieta = 0.;
	std::vector<double>* TrkSeta = 0.;
	std::vector<double>* TrkKphi = 0.;
	std::vector<double>* Trkpiphi = 0.;
	std::vector<double>* TrkSphi = 0.;
	std::vector<double>* TrkScharge = 0.;
	std::vector<double>* D0fromDSsXY_vec = 0.;
	std::vector<double>* D0mass = 0.;
	std::vector<double>* D0mass = 0.;
	std::vector<double>* D0mass = 0.;

	
	std::vector<double> vectorInvariantMass_Dstar = 0.;
	std::vector<double> vectorInvariantMass_D0 = 0.;

	//counters for file f1 (pythia)
	int countTotalEvents = 0;

	int countTotalTracks = 0;
	int countTracksAfterTrigger = 0;
	int countTracksHasTrackDetails = 0;
	int countTracksChargeZero = 0;
	int countTracksHighPurity = 0;
	int countTracksPDG211 = 0;
	int countTracksPtZeroFive = 0;
	int countTracksNumberOfHits2 = 0;
	int countTracksDxyThree = 0;
	int countTracksDzThree = 0;
	int countTrackSlowPionCandidates = 0;
	int countObservation = 0;
	int countTracksPtZeroSix = 0;
	int countTracksChiTwoFive = 0;
	int countTracksNumberOfHits5 = 0;
	int countTracksNumberOfPixelHits2 = 0;
	int countTracksDzOne = 0;
	int countTracksDxyZeroOne = 0;
	int countTrackKaonPionCandidates = 0;

	int countD0AfterLorentzVector = 0;
	int countDsMinusD0Zerothree = 0;
	int countTransientTrackOfpiK = 0;
	int countPointingcosPhi = 0;
	int countSignificance = 0;
	int countD0pTThree = 0;
	int countD0Candidates = 0;
	int countDsAfterLorentzVector = 0;
	int countDsMinusD0 = 0;
	int countDsCandidates = 0;

	//double M = 0.;
	//double Pt = 0.;
	//double Eta = 0.;
	//double Rapidity = 0.;
	
	//*****Creating Histgrams******************************************************	

	//Histogramas cinematics quantities of the D0
	TH1D *D0mass_Histo = new TH1D("D0mass_Histo","D0mass_Histo",100,1.7,2.0);
	D0mass_Histo->SetTitle("Invariant Mass of the D0 ; Mass [GeV] ; Events ");
	D0mass_Histo->SetName("D0mass_Histo");

	TH1F *D0eta_Histo = new TH1F("D0eta_Histo","D0eta_Histo",100,-3,3);
	D0eta_Histo->SetTitle("Pseudo-rapidity distribuition of the D0 ; #eta ; Events ");
	D0eta_Histo->SetName("D0eta_Histo");

	TH1F *D0phi_Histo = new TH1F("D0phi_Histo","D0phi_Histo",100,-5,5);
	D0phi_Histo->SetTitle("#Phi distribuition of the D0 ; #Phi ; Events ");
	D0phi_Histo->SetName("D0phi_Histo");

	TH1F *D0pt_Histo = new TH1F("D0pt_Histo","D0pt_Histo",100,0,40);
	D0pt_Histo->SetTitle("pT distribuition of the D0 ; pT [GeV] ; Events ");
	D0pt_Histo->SetName("D0pt_Histo");

	//Histogramas cinematics quantities of the D*
	TH1D *Dsmass_Histo = new TH1D("Dsmass_Histo","Dsmass_Histo",100,1.85,2.15);
	Dsmass_Histo->SetTitle("Invariant Mass of the D* ; Mass [GeV] ; Events ");
	Dsmass_Histo->SetName("Dsmass_Histo");

	TH1F *Dseta_Histo = new TH1F("Dseta_Histo","Dseta_Histo",100,-3,3);
	Dseta_Histo->SetTitle("Pseudo-rapidity distribuition of the D* ; #eta ; Events ");
	Dseta_Histo->SetName("Dseta_Histo");

	TH1F *Dsphi_Histo = new TH1F("Dsphi_Histo","Dsphi_Histo",100,-5,5);
	Dsphi_Histo->SetTitle("#Phi distribuition of the D* ; #Phi ; Events ");
	Dsphi_Histo->SetName("Dsphi_Histo");

	TH1F *Dspt_Histo = new TH1F("Dspt_Histo","DspT_Histo",100,0,40);
	Dspt_Histo->SetTitle("pT distribuition of the D* ; pT [GeV] ; Events ");
	Dspt_Histo->SetName("Dspt_Histo");

	//Histogramas cinematics quantities of the SlowPion
	TH1F *TrkSpt_Histo = new TH1F("TrkSpt_Histo","TrkSpt_Histo",100,0,2.5);
	TrkSpt_Histo->SetTitle("pt distribuition of the SlowPion ; pT [GeV] ; Events ");
	TrkSpt_Histo->SetName("TrkSpt_Histo");

	TH1F *TrkSnhits_Histo = new TH1F("TrkSnhits_Histo","TrkSnhits_Histo",100,0,45);
	TrkSnhits_Histo->SetTitle("nhits distribuition of the SlowPion ; Number of hits ; Events ");
	TrkSnhits_Histo->SetName("TrkSnhits_Histo");

	TH1F *TrkSdxy_Histo = new TH1F("TrkSdxy_Histo","TrkSdxy_Histo",100,-2,2);
	TrkSdxy_Histo->SetTitle("dxy distribuition of the SlowPion ; dx [cm] ; Events ");
	TrkSdxy_Histo->SetName("TrkSdxy_Histo");

	TH1F *TrkSdz_Histo = new TH1F("TrkSdz_Histo","TrkSdz_Histo",100,-4,4);
	TrkSdz_Histo->SetTitle("dz distribuition of the SlowPion ; dz [cm] ; Events ");
	TrkSdz_Histo->SetName("TrkSdz_Histo");

	TH1F *TrkSeta_Histo = new TH1F("TrkSeta_Histo","TrkSeta_Histo",100,-2.5,2.5);
	TrkSeta_Histo->SetTitle("eta distribuition of the SlowPion ; #eta ; Events ");
	TrkSeta_Histo->SetName("TrkSeta_Histo");

	TH1F *TrkSphi_Histo = new TH1F("TrkSphi_Histo","TrkSphi_Histo",100,-4,4);
	TrkSphi_Histo->SetTitle("Phi distribuition of the SlowPion ; #Phi ; Events ");
	TrkSphi_Histo->SetName("TrkSphi_Histo");

	TH1F *TrkSchi2_Histo = new TH1F("TrkSchi2_Histo","TrkSchi2_Histo",10,0,3);
	TrkSchi2_Histo->SetTitle("Chi2 distribuition of the SlowPion ; #Chi^{2} ; Events ");
	TrkSchi2_Histo->SetName("TrkSchi2_Histo");

	//Histogramas cinematics quantities of the Pion
	TH1F *Trkpipt_Histo = new TH1F("Trkpipt_Histo","Trkpipt_Histo",100,0,10);
	Trkpipt_Histo->SetTitle("pt distribuition of the Pion ; pT [GeV] ; Events ");
	Trkpipt_Histo->SetName("Trkpipt_Histo");

	TH1F *Trkpinhits_Histo = new TH1F("Trkpinhits_Histo","Trkpinhits_Histo",100,0,40);
	Trkpinhits_Histo->SetTitle("nhits distribuition of the Pion ; Number of hits ; Events ");
	Trkpinhits_Histo->SetName("Trkpinhits_Histo");

	TH1F *Trkpidxy_Histo = new TH1F("Trkpidxy_Histo","Trkpidxy_Histo",100,-0.8,0.8);
	Trkpidxy_Histo->SetTitle("dxy distribuition of the Pion ; dx [cm] ; Events ");
	Trkpidxy_Histo->SetName("Trkpidxy_Histo");

	TH1F *Trkpidz_Histo = new TH1F("Trkpidz_Histo","Trkpidz_Histo",100,-1.5,1.5);
	Trkpidz_Histo->SetTitle("dz distribuition of the Pion ; dz [cm] ; Events ");
	Trkpidz_Histo->SetName("Trkpidz_Histo");

	TH1F *Trkpieta_Histo = new TH1F("Trkpieta_Histo","Trkpieta_Histo",100,-2.5,2.5);
	Trkpieta_Histo->SetTitle("eta distribuition of the Pion ; #eta ; Events ");
	Trkpieta_Histo->SetName("Trkpieta_Histo");

	TH1F *Trkpiphi_Histo = new TH1F("Trkpiphi_Histo","Trkpiphi_Histo",100,-4,4);
	Trkpiphi_Histo->SetTitle("phi distribuition of the Pion ; #Phi ; Events ");
	Trkpiphi_Histo->SetName("Trkpiphi_Histo");

	TH1F *Trkpichi2_Histo = new TH1F("Trkpichi2_Histo","Trkpichi2_Histo",10,0,3);
	Trkpichi2_Histo->SetTitle("Chi2 distribuition of the Pion ; #Chi^{2} ; Events ");
	Trkpichi2_Histo->SetName("Trkpichi2_Histo");

	//Histogramas cinematics quantities of the Kaon
	TH1F *TrkKpt_Histo = new TH1F("TrkKpt_Histo","TrkKpt_Histo",100,0,10);
	TrkKpt_Histo->SetTitle("pt distribuition of the Kaon ; pT [GeV] ; Events ");
	TrkKpt_Histo->SetName("TrkKpt_Histo");

	TH1F *TrkKnhits_Histo = new TH1F("TrkKnhits_Histo","TrkKnhits_Histo",100,0,40);
	TrkKnhits_Histo->SetTitle("nhits distribuition of the Kaon ; Number of hits ; Events ");
	TrkKnhits_Histo->SetName("TrkKnhits_Histo");

	TH1F *TrkKdxy_Histo = new TH1F("TrkKdxy_Histo","TrkKdxy_Histo",100,-0.8,0.8);
	TrkKdxy_Histo->SetTitle("dxy distribuition of the Kaon ; dx [cm] ; Events ");
	TrkKdxy_Histo->SetName("TrkKdxy_Histo");

	TH1F *TrkKdz_Histo = new TH1F("TrkKdz_Histo","TrkKdz_Histo",100,-1.5,1.5);
	TrkKdz_Histo->SetTitle("dz distribuition of the Kaon ; dz [cm] ; Events ");
	TrkKdz_Histo->SetName("TrkKdz_Histo");

	TH1F *TrkKeta_Histo = new TH1F("TrkKeta_Histo","TrkKeta_Histo",100,-2.5,2.5);
	TrkKeta_Histo->SetTitle("eta distribuition of the Kaon ; #eta ; Events ");
	TrkKeta_Histo->SetName("TrkKeta_Histo");

	TH1F *TrkKphi_Histo = new TH1F("TrkKphi_Histo","TrkKphi_Histo",100,-4,4);
	TrkKphi_Histo->SetTitle("phi distribuition of the Kaon ; #Phi ; Events ");
	TrkKphi_Histo->SetName("TrkKphi_Histo");

	TH1F *TrkKchi2_Histo = new TH1F("TrkKchi2_Histo","TrkKchi2_Histo",10,0,3);
	TrkKchi2_Histo->SetTitle("Chi2 distribuition of the Kaon ; #Chi^{2} ; Events ");
	TrkKchi2_Histo->SetName("TrkKchi2_Histo");

	//End Histograms
	//-----------------------------------------------------------------------------

	//-------Reading the root file and the tree-------------------------------------	
	TFile *f1 = new TFile("D0DstarData_parcialfiles.root");
	TTree *t1 = (TTree*)f1->Get("analysis/data");

	//---------------------------------------------------------------------------------
	// addressing the memory to vector and variables for file

	//For Variables
	//TBranch *b_Total_Events = t1->GetBranch("Total_Events");
	//b_Total_Events->SetAddress(&Total_Events);

	if(Conters)
	{
	//counters
	TBranch *b_Total_Events = t1->GetBranch("Total_Events");
	b_Total_Events->SetAddress(&Total_Events);
	TBranch *b_TotalTracks = t1->GetBranch("TotalTracks");
	b_TotalTracks->SetAddress(&TotalTracks);
	TBranch *b_TracksAfterTrigger = t1->GetBranch("TracksAfterTrigger");
	b_TracksAfterTrigger->SetAddress(&TracksAfterTrigger);
	TBranch *b_TracksHasTrackDetails = t1->GetBranch("TracksHasTrackDetails");
	b_TracksHasTrackDetails->SetAddress(&TracksHasTrackDetails);
	TBranch *b_TracksChargeZero = t1->GetBranch("TracksChargeZero");
	b_TracksChargeZero->SetAddress(&TracksChargeZero);
	TBranch *b_TracksHighPurity = t1->GetBranch("TracksHighPurity");
	b_TracksHighPurity->SetAddress(&TracksHighPurity);
	TBranch *b_TracksPDG211 = t1->GetBranch("TracksPDG211");
	b_TracksPDG211->SetAddress(&TracksPDG211);
	TBranch *b_TracksPtZeroFive = t1->GetBranch("TracksPtZeroFive");
	b_TracksPtZeroFive->SetAddress(&TracksPtZeroFive);
	TBranch *b_TracksNumberOfHits2 = t1->GetBranch("TracksNumberOfHits2");
	b_TracksNumberOfHits2->SetAddress(&TracksNumberOfHits2);
	TBranch *b_TracksDxyThree = t1->GetBranch("TracksDxyThree");
	b_TracksDxyThree->SetAddress(&TracksDxyThree);
	TBranch *b_TracksDzThree = t1->GetBranch("TracksDzThree");
	b_TracksDzThree->SetAddress(&TracksDzThree);
	TBranch *b_TrackSlowPionCandidates = t1->GetBranch("TrackSlowPionCandidates");
	b_TrackSlowPionCandidates->SetAddress(&TrackSlowPionCandidates);
	TBranch *b_Observation = t1->GetBranch("Observation");
	b_Observation->SetAddress(&Observation);
	TBranch *b_TracksPtZeroSix = t1->GetBranch("TracksPtZeroSix");
	b_TracksPtZeroSix->SetAddress(&TracksPtZeroSix);
	TBranch *b_TracksChiTwoFive = t1->GetBranch("TracksChiTwoFive");
	b_TracksChiTwoFive->SetAddress(&TracksChiTwoFive);
	TBranch *b_TracksNumberOfHits5 = t1->GetBranch("TracksNumberOfHits5");
	b_TracksNumberOfHits5->SetAddress(&TracksNumberOfHits5);
	TBranch *b_TracksNumberOfPixelHits2 = t1->GetBranch("TracksNumberOfPixelHits2");
	b_TracksNumberOfPixelHits2->SetAddress(&TracksNumberOfPixelHits2);
	TBranch *b_TracksDzOne = t1->GetBranch("TracksDzOne");
	b_TracksDzOne->SetAddress(&TracksDzOne);
	TBranch *b_TracksDxyZeroOne = t1->GetBranch("TracksDxyZeroOne");
	b_TracksDxyZeroOne->SetAddress(&TracksDxyZeroOne);
	TBranch *b_TrackKaonPionCandidates = t1->GetBranch("TrackKaonPionCandidates");
	b_TrackKaonPionCandidates->SetAddress(&TrackKaonPionCandidates);

	TBranch *b_D0AfterLorentzVector = t1->GetBranch("D0AfterLorentzVector");
	b_D0AfterLorentzVector->SetAddress(&D0AfterLorentzVector);
	TBranch *b_DsMinusD0Zerothree = t1->GetBranch("DsMinusD0Zerothree");
	b_DsMinusD0Zerothree->SetAddress(&DsMinusD0Zerothree); 
	TBranch *b_TransientTrackOfpiK = t1->GetBranch("TransientTrackOfpiK");
	b_TransientTrackOfpiK->SetAddress(&TransientTrackOfpiK);
	TBranch *b_PointingcosPhi = t1->GetBranch("PointingcosPhi");
	b_PointingcosPhi->SetAddress(&PointingcosPhi);
	TBranch *b_Significance = t1->GetBranch("Significance");
	b_Significance->SetAddress(&Significance);
	TBranch *b_D0pTThree = t1->GetBranch("D0pTThree");
	b_D0pTThree->SetAddress(&D0pTThree);
	TBranch *b_D0Candidates = t1->GetBranch("D0Candidates");
	b_D0Candidates->SetAddress(&D0Candidates);
	TBranch *b_DsAfterLorentzVector = t1->GetBranch("DsAfterLorentzVector");
	b_DsAfterLorentzVector->SetAddress(&DsAfterLorentzVector);
	TBranch *b_DsMinusD0 = t1->GetBranch("DsMinusD0");
	b_DsMinusD0->SetAddress(&DsMinusD0);
	TBranch *b_DsCandidates = t1->GetBranch("DsCandidates");
	b_DsCandidates->SetAddress(&DsCandidates);
	}

	if(VectorsHistograms)
	{		
	//For Vectors	
	TBranch *b_D0mass = t1->GetBranch("D0mass");
	b_D0mass->SetAddress(&D0mass);
	TBranch *b_D0eta = t1->GetBranch("D0eta");
	b_D0eta->SetAddress(&D0eta);
	TBranch *b_D0phi = t1->GetBranch("D0phi");
	b_D0phi->SetAddress(&D0phi);
	TBranch *b_D0pt = t1->GetBranch("D0pt");
	b_D0pt->SetAddress(&D0pt);  

	TBranch *b_Dsmass = t1->GetBranch("Dsmass");
	b_Dsmass->SetAddress(&Dsmass);
	TBranch *b_Dseta = t1->GetBranch("Dseta");
	b_Dseta->SetAddress(&Dseta);
	TBranch *b_Dsphi = t1->GetBranch("Dsphi");
	b_Dsphi->SetAddress(&Dsphi);
	TBranch *b_Dspt = t1->GetBranch("Dspt");
	b_Dspt->SetAddress(&Dspt);

	TBranch *b_TrkSpt = t1->GetBranch("TrkSpt");
	b_TrkSpt->SetAddress(&TrkSpt);
	TBranch *b_TrkSnhits = t1->GetBranch("TrkSnhits");
	b_TrkSnhits->SetAddress(&TrkSnhits);
	TBranch *b_TrkSdxy = t1->GetBranch("TrkSdxy");
	b_TrkSdxy->SetAddress(&TrkSdxy);
	TBranch *b_TrkSdz = t1->GetBranch("TrkSdz");
	b_TrkSdz->SetAddress(&TrkSdz);
	TBranch *b_TrkSeta = t1->GetBranch("TrkSeta");
	b_TrkSeta->SetAddress(&TrkSeta);
	TBranch *b_TrkSphi = t1->GetBranch("TrkSphi");
	b_TrkSphi->SetAddress(&TrkSphi);
	TBranch *b_TrkSchi2 = t1->GetBranch("TrkSchi2");
	b_TrkSchi2->SetAddress(&TrkSchi2);

	TBranch *b_Trkpipt = t1->GetBranch("Trkpipt");
	b_Trkpipt->SetAddress(&Trkpipt);
	TBranch *b_Trkpinhits = t1->GetBranch("Trkpinhits");
	b_Trkpinhits->SetAddress(&Trkpinhits);
	TBranch *b_Trkpidxy = t1->GetBranch("Trkpidxy");
	b_Trkpidxy->SetAddress(&Trkpidxy);
	TBranch *b_Trkpidz = t1->GetBranch("Trkpidz");
	b_Trkpidz->SetAddress(&Trkpidz);
	TBranch *b_Trkpieta = t1->GetBranch("Trkpieta");
	b_Trkpieta->SetAddress(&Trkpieta);
	TBranch *b_Trkpiphi = t1->GetBranch("Trkpiphi");
	b_Trkpiphi->SetAddress(&Trkpiphi);
	TBranch *b_Trkpichi2 = t1->GetBranch("Trkpichi2");
	b_Trkpichi2->SetAddress(&Trkpichi2);

	TBranch *b_TrkKpt = t1->GetBranch("TrkKpt");
	b_TrkKpt->SetAddress(&TrkKpt);
	TBranch *b_TrkKnhits = t1->GetBranch("TrkKnhits");
	b_TrkKnhits->SetAddress(&TrkKnhits);
	TBranch *b_TrkKdxy = t1->GetBranch("TrkKdxy");
	b_TrkKdxy->SetAddress(&TrkKdxy);
	TBranch *b_TrkKdz = t1->GetBranch("TrkKdz");
	b_TrkKdz->SetAddress(&TrkKdz);
	TBranch *b_TrkKeta = t1->GetBranch("TrkKeta");
	b_TrkKeta->SetAddress(&TrkKeta);
	TBranch *b_TrkKphi = t1->GetBranch("TrkKphi");
	b_TrkKphi->SetAddress(&TrkKphi);
	TBranch *b_TrkKchi2 = t1->GetBranch("TrkKchi2");
	b_TrkKchi2->SetAddress(&TrkKchi2);
	}

	//**********************************************************		
	//Reading Number of tree entries for file f1
	Long64_t nentries = t1->GetEntries();
	cout<< "Number of tree entries: "<< nentries <<std::endl;

	Long64_t GetEntriesFast = t1->GetEntriesFast();
	//cout<< "GetEntriesFast: "<< GetEntriesFast <<std::endl;

	Long64_t nbytes = 0, nb = 0, i=0;
	for (Long64_t jentry=0; jentry < 2000000; jentry++) //loop tree entries for file f1
	{
      	Long64_t ientry = t1->LoadTree(jentry);
      	//std::cout << "nentries " << nentries << " ientry " << ientry << " jentry " << jentry <<std::endl;
      	if (ientry < 0) break;

		//if ((jentry/nentries) % (nentries*0.01) == 0) cout <<"*****"<< (jentry*100)/nentries << "per cent done***" << endl;
		//{
			
		//}

		double percent = (jentry*100)/2000000;
		if (jentry % 10000 == 0) cout <<"*****"<< percent << "per cent done***" << endl;
		//cout <<"*****entrada"<< jentry << "***" << endl;

		if(Conters)
		{
		//counters
		b_Total_Events->GetEntry(ientry);				
		b_TotalTracks->GetEntry(ientry);
		b_TracksAfterTrigger->GetEntry(ientry);
		b_TracksHasTrackDetails->GetEntry(ientry);
		b_TracksChargeZero->GetEntry(ientry);
		b_TracksHighPurity->GetEntry(ientry);
		b_TracksPDG211->GetEntry(ientry);
		b_TracksPtZeroFive->GetEntry(ientry);
		b_TracksNumberOfHits2->GetEntry(ientry);
		b_TracksDxyThree->GetEntry(ientry);
		b_TracksDzThree->GetEntry(ientry);
		b_TrackSlowPionCandidates->GetEntry(ientry);
		b_TracksPtZeroSix->GetEntry(ientry); 
		b_Observation->GetEntry(ientry);
		b_TracksChiTwoFive->GetEntry(ientry);
		b_TracksNumberOfHits5->GetEntry(ientry);
		b_TracksNumberOfPixelHits2->GetEntry(ientry);
		b_TracksDxyZeroOne->GetEntry(ientry);
		b_TracksDzOne->GetEntry(ientry);
		b_TrackKaonPionCandidates->GetEntry(ientry);

		b_D0AfterLorentzVector->GetEntry(ientry);
		b_DsMinusD0Zerothree->GetEntry(ientry);
		b_TransientTrackOfpiK->GetEntry(ientry);
		b_PointingcosPhi->GetEntry(ientry);
		b_Significance->GetEntry(ientry);
		b_D0pTThree->GetEntry(ientry);
		b_D0Candidates->GetEntry(ientry);
		b_DsMinusD0->GetEntry(ientry);
		b_DsCandidates->GetEntry(ientry);
		b_DsAfterLorentzVector->GetEntry(ientry);

		countTotalEvents += Total_Events;
		countTotalTracks += TotalTracks;
		countTracksAfterTrigger += TracksAfterTrigger;
		countTracksHasTrackDetails += TracksHasTrackDetails;
		countTracksChargeZero += TracksChargeZero;
		countTracksHighPurity += TracksHighPurity;
		countTracksPDG211 += TracksPDG211;
		countTracksPtZeroFive += TracksPtZeroFive;
		countTracksNumberOfHits2 += TracksNumberOfHits2;
		countTracksDxyThree += TracksDxyThree;
		countTracksDzThree += TracksDzThree;
		countTrackSlowPionCandidates += TrackSlowPionCandidates;
		countTracksPtZeroSix += TracksPtZeroSix;
		countObservation += Observation;
		countTracksChiTwoFive += TracksChiTwoFive;
		countTracksNumberOfHits5 += TracksNumberOfHits5;
		countTracksNumberOfPixelHits2 += TracksNumberOfPixelHits2;
		countTracksDxyZeroOne += TracksDxyZeroOne;
		countTracksDzOne += TracksDzOne;
		countTrackKaonPionCandidates += TrackKaonPionCandidates;

		countD0AfterLorentzVector += D0AfterLorentzVector;
		countDsMinusD0Zerothree += DsMinusD0Zerothree;
		countTransientTrackOfpiK += TransientTrackOfpiK;
		countPointingcosPhi += PointingcosPhi;
		countSignificance += Significance;
		countD0pTThree += D0pTThree;
		countD0Candidates += D0Candidates;
		countDsAfterLorentzVector += DsAfterLorentzVector;
		countDsMinusD0 += DsMinusD0;
		countDsCandidates += DsCandidates;
		}//End Flag Conters		

		if(VectorsHistograms)
		{
		//vectors
       	b_D0mass->GetEntry(ientry);
		b_D0eta->GetEntry(ientry);
		b_D0phi->GetEntry(ientry);
		b_D0pt->GetEntry(ientry);

		b_Dsmass->GetEntry(ientry);
		b_Dseta->GetEntry(ientry);	
		b_Dsphi->GetEntry(ientry); 
		b_Dspt->GetEntry(ientry);

		if(Meson)
		{
		b_TrkSpt->GetEntry(ientry);
		b_TrkSnhits->GetEntry(ientry);
		b_TrkSdxy->GetEntry(ientry);
		b_TrkSdz->GetEntry(ientry);
		b_TrkSeta->GetEntry(ientry);
		b_TrkSphi->GetEntry(ientry);
		b_TrkSchi2->GetEntry(ientry);

		b_Trkpipt->GetEntry(ientry);
		b_Trkpinhits->GetEntry(ientry);
		b_Trkpidxy->GetEntry(ientry);
		b_Trkpidz->GetEntry(ientry);
		b_Trkpieta->GetEntry(ientry);
		b_Trkpiphi->GetEntry(ientry);
		b_Trkpichi2->GetEntry(ientry);

		b_TrkKpt->GetEntry(ientry);
		b_TrkKnhits->GetEntry(ientry);
		b_TrkKdxy->GetEntry(ientry);
		b_TrkKdz->GetEntry(ientry);
		b_TrkKeta->GetEntry(ientry);
		b_TrkKphi->GetEntry(ientry);
		b_TrkKchi2->GetEntry(ientry);
		}

		//For D0
		for(Long64_t i=0; i < D0mass->size(); i++)
		{  
			D0mass_Histo->Fill(D0mass->at(i));
			vectorInvariantMass_D0.push_back(D0mass->at(i));
  		}
		for(Long64_t i=0; i < D0eta->size(); i++)
		{  
			D0eta_Histo->Fill(D0eta->at(i));
  		}
		for(Long64_t i=0; i < D0phi->size(); i++)
		{  
			D0phi_Histo->Fill(D0phi->at(i));
  		}
		for(Long64_t i=0; i < D0phi->size(); i++)
		{  
			D0pt_Histo->Fill(D0pt->at(i));
  		}

		//For D*
		for(Long64_t i=0; i < Dsmass->size(); i++)
		{  
			Dsmass_Histo->Fill(Dsmass->at(i));
			vectorInvariantMass_Dstar.push_back(Dsmass->at(i));
  		}
		for(Long64_t i=0; i < Dseta->size(); i++)
		{  
			Dseta_Histo->Fill(Dseta->at(i));
  		}
		for(Long64_t i=0; i < Dsphi->size(); i++)
		{  
			Dsphi_Histo->Fill(Dsphi->at(i));
  		}
		for(Long64_t i=0; i < Dsphi->size(); i++)
		{  
			Dspt_Histo->Fill(Dspt->at(i));
  		}

		if(Meson)
		{
		//For SlowPion
		for(Long64_t i=0; i < TrkSpt->size(); i++)
		{  
			TrkSpt_Histo->Fill(TrkSpt->at(i));
  		}

		for(Long64_t i=0; i < TrkSnhits->size(); i++)
		{  
			TrkSnhits_Histo->Fill(TrkSnhits->at(i));
  		}

		for(Long64_t i=0; i < TrkSdxy->size(); i++)
		{  
			TrkSdxy_Histo->Fill(TrkSdxy->at(i));
			//std::cout << "TrkSdxy->at(i): " << TrkSdxy->at(i) << endl;
  		}

		for(Long64_t i=0; i < TrkSdz->size(); i++)
		{  
			TrkSdz_Histo->Fill(TrkSdz->at(i));
  		}

		for(Long64_t i=0; i < TrkSeta->size(); i++)
		{  
			TrkSeta_Histo->Fill(TrkSeta->at(i));
  		}

		for(Long64_t i=0; i < TrkSphi->size(); i++)
		{  
			TrkSphi_Histo->Fill(Dsphi->at(i));
  		}
	
		for(Long64_t i=0; i < TrkSchi2->size(); i++)
		{  
			TrkSchi2_Histo->Fill(TrkSchi2->at(i));
  		}

		//For Pion
		for(Long64_t i=0; i < Trkpipt->size(); i++)
		{  
			Trkpipt_Histo->Fill(Trkpipt->at(i));
  		}

		for(Long64_t i=0; i < Trkpinhits->size(); i++)
		{  
			Trkpinhits_Histo->Fill(Trkpinhits->at(i));
  		}

		for(Long64_t i=0; i < Trkpidxy->size(); i++)
		{  
			Trkpidxy_Histo->Fill(Trkpidxy->at(i));
  		}

		for(Long64_t i=0; i < Trkpidz->size(); i++)
		{  
			Trkpidz_Histo->Fill(Trkpidz->at(i));
  		}

		for(Long64_t i=0; i < Trkpieta->size(); i++)
		{  
			Trkpieta_Histo->Fill(Trkpieta->at(i));
  		}

		for(Long64_t i=0; i < Trkpiphi->size(); i++)
		{  
			Trkpiphi_Histo->Fill(Trkpiphi->at(i));
  		}

		for(Long64_t i=0; i < Trkpichi2->size(); i++)
		{  
			Trkpichi2_Histo->Fill(Trkpichi2->at(i));
  		}

		//For Kon
		for(Long64_t i=0; i < TrkKpt->size(); i++)
		{  
			TrkKpt_Histo->Fill(TrkKpt->at(i));
  		}

		for(Long64_t i=0; i < TrkKnhits->size(); i++)
		{  
			TrkKnhits_Histo->Fill(TrkKnhits->at(i));
  		}

		for(Long64_t i=0; i < TrkKdxy->size(); i++)
		{  
			TrkKdxy_Histo->Fill(TrkKdxy->at(i));
  		}

		for(Long64_t i=0; i < TrkKdz->size(); i++)
		{  
			TrkKdz_Histo->Fill(TrkKdz->at(i));
  		}

		for(Long64_t i=0; i < TrkKeta->size(); i++)
		{  
			TrkKeta_Histo->Fill(TrkKeta->at(i));
  		}

		for(Long64_t i=0; i < TrkKphi->size(); i++)
		{  
			TrkKphi_Histo->Fill(TrkKphi->at(i));
  		}

		for(Long64_t i=0; i < TrkKchi2->size(); i++)
		{  
			TrkKchi2_Histo->Fill(TrkKchi2->at(i));
  		}
		}//End Flag Meson
		}//End Flag VectorsHistograms
			
	}//End loop tree entries for file f1

	TFile f_analysis("DMesonsMass.root","recreate"); //Creates root file
	TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	t_analysis.Branch("vectorInvariantMass_Dstar",&vectorInvariantMass_Dstar); //Creates a branch
	t_analysis.Branch("vectorInvariantMass_D0",&vectorInvariantMass_D0); //Creates a branch
	t_analysis.Fill();
	Dsmass_Histo->Write(); //Write() if a file is open, this function writes a root objectics on it.
	D0mass_Histo->Write(); //Write() if a file is open, this function writes a root objectics on it.
	t_analysis.Write();  //Write in the root file

	//TFile f_analysis("Jpsi.root","recreate"); //Creates root file
	//TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	//t_analysis.Branch("vectorInvariantMassJpsi",&vectorInvariantMassJpsi); //Creates a branch
	//t_analysis.Fill();
	//h_DimuonsOppositeChargeEtaJpsi_M2_Double.Write(); //Write() if a file is open, this function writes a root objectics on it.
	//t_analysis.Write();  //Write in the root file


	cout << "countTotalEvents: "<< countTotalEvents << endl;
	cout << "        " << endl;	
	cout << "=============Tracks=========================================== " << endl;
	cout << "countTotalTracks: "<< countTotalTracks << endl;
	cout << "countTracksAfterTrigger: "<< countTracksAfterTrigger << endl;
	cout << "countTracksHasTrackDetails: "<< countTracksHasTrackDetails << endl;
	cout << "countTracksChargeZero: "<< countTracksChargeZero << endl;
	cout << "countTracksHighPurity: "<< countTracksHighPurity << endl;
	cout << "countTracksPDG211: "<< countTracksPDG211 << endl;
	cout << "=============Criteria to SlowPion Tracks=======================" << endl;
	cout << "countTracksPtZeroFive: "<< countTracksPtZeroFive << endl;
	cout << "countTracksNumberOfHits2: "<< countTracksNumberOfHits2 << endl;
	cout << "countTracksDxyThree: "<< countTracksDxyThree << endl;
	cout << "countTracksDzThree: "<< countTracksDzThree << endl;
	cout << "countTrackSlowPionCandidates: "<< countTrackSlowPionCandidates << endl;
	cout << "Observation: "<< Observation << endl;
	cout << "=============Criteria to Pions and Kaons Tracks=======================" << endl;
	cout << "countTracksPtZeroSix: "<< countTracksPtZeroSix << endl;
	cout << "countTracksChiTwoFive: "<< countTracksChiTwoFive << endl;
	cout << "countTracksNumberOfHits5: "<< countTracksNumberOfHits5 << endl;
	cout << "countTracksNumberOfPixelHits2: "<< countTracksNumberOfPixelHits2 << endl;
	cout << "countTracksDxyZeroOne: "<< countTracksDxyZeroOne << endl;
	cout << "countTracksDzOne: "<< countTracksDzOne << endl;
	cout << "countTrackKaonPionCandidates: "<< countTrackKaonPionCandidates << endl;
	cout << "=============Criteria to D0 and D* =======================" << endl;
	cout << "countD0AfterLorentzVector: "<< countD0AfterLorentzVector << endl;
	cout << "countDsMinusD0Zerothree: "<< countDsMinusD0Zerothree << endl;
	cout << "countTransientTrackOfpiK: "<< countTransientTrackOfpiK << endl;
	cout << "countPointingcosPhi: "<< countPointingcosPhi << endl;
	cout << "countSignificance: "<< countSignificance << endl;
	cout << "countD0pTThree: "<< countD0pTThree << endl;
	cout << "countD0Candidates: "<< countD0Candidates << endl;
	cout << "countDsAfterLorentzVector: "<< countDsAfterLorentzVector << endl;
	cout << "countDsMinusD0: "<< countDsMinusD0 << endl;
	cout << "countDsCandidates: "<< countDsCandidates << endl;
	
	if(VectorsHistograms)
	{
    //========Creating Canvas=================================================================	
	TCanvas* c0 = new TCanvas("c0","Canvas 0 - behavior of the D",1200,600);
	c0->Divide(2);
	c0->cd(1);
	//D0mass_Histo->SetLineColor(kBlue);
	D0mass_Histo->SetMarkerStyle(7);
	//D0mass_Histo->SetMarkerStyle(21);
	//D0mass_Histo->SetStats(0);
	//D0mass_Histo->SetMarkerColor(kBlue);
	//D0mass_Histo->SetFillColor(kBlue);
	D0mass_Histo->Sumw2();
	D0mass_Histo->Draw("e1p");
	//D0mass_Histo->Draw("e0");
	//D0mass_Histo->Draw("e1pSAME");
	/*TLegend* leg_D0mass = new TLegend(0.5,0.5,0.75,0.65);
   	leg_D0mass->SetFillColor(kWhite);
	leg_D0mass->SetFillStyle(1001);
	leg_D0mass->SetBorderSize(0);
	leg_D0mass->AddEntry(D0mass_Histo,"pT > 3 GeV","e1pSAME");
	leg_D0mass->AddEntry(D0mass_Histo,"D0mass","L");
	leg_D0mass->Draw();*/
	//-------------------------------------------------------------------------
	c0->cd(2);
	Dsmass_Histo->SetLineColor(kBlue);
	Dsmass_Histo->SetMarkerStyle(7);
    Dsmass_Histo->Sumw2();
	Dsmass_Histo->Draw("e1p");
	//Dsmass_Histo->Draw("e1pSAME");
	/*TLegend* leg_Dsmass = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Dsmass->SetFillColor(kWhite);
	leg_Dsmass->SetFillStyle(1001);
	leg_Dsmass->SetBorderSize(0);
	//leg_Dsmass->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Dsmass->AddEntry(Dsmass_Histo,"Dsmass","L");
	leg_Dsmass->Draw();*/
	c0->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas0.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c1 = new TCanvas("c1","Canvas 1 - behavior of the D",1200,600);
	c1->Divide(2);
	c1->cd(1);
	//D0eta_Histo->SetLineColor(kBlue);
	D0pt_Histo->SetMarkerStyle(7);
	D0pt_Histo->Draw("HIST");
	//D0pt_Histo->Draw("e1pSAME");
	TLegend* leg_D0pt = new TLegend(0.5,0.5,0.65,0.6);
   	leg_D0pt->SetFillColor(kWhite);
	leg_D0pt->SetFillStyle(1001);
	leg_D0pt->SetBorderSize(0);
	leg_D0pt->AddEntry(D0pt_Histo,"pT > 3 GeV","L" );
	//leg_D0pt->AddEntry(D0eta_Histo,"D0eta","L");
	leg_D0pt->Draw();
	//-------------------------------------------------------------------------
	c1->cd(2);
	//Dseta_Histo->SetLineColor(kBlue);
	Dspt_Histo->SetMarkerStyle(7);
	Dspt_Histo->Draw("HIST");
	//Dspt_Histo->Draw("e1pSAME");	
	/*TLegend* leg_Dspt = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Dspt->SetFillColor(kWhite);
	leg_Dspt->SetFillStyle(1001);
	leg_Dspt->SetBorderSize(0);
	//leg_Dspt->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Dspt->AddEntry(Dspt_Histo,"Dseta","L");
	leg_Dspt->Draw();*/
	c1->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas1.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c2 = new TCanvas("c2","Canvas 2 - behavior of the D",1200,600);
	c2->Divide(2);
	c2->cd(1);
	//D0eta_Histo->SetLineColor(kBlue);
	D0eta_Histo->SetMarkerStyle(7);
	D0eta_Histo->Draw("HIST");
	//D0eta_Histo->Draw("e1pSAME");
	/*TLegend* leg_D0eta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_D0eta->SetFillColor(kWhite);
	leg_D0eta->SetFillStyle(1001);
	leg_D0eta->SetBorderSize(0);
	//leg_D0eta->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_D0eta->AddEntry(D0eta_Histo,"D0eta","L");
	leg_D0eta->Draw();*/
	//-------------------------------------------------------------------------
	c2->cd(2);
	//Dseta_Histo->SetLineColor(kBlue);
	Dseta_Histo->SetMarkerStyle(7);
	Dseta_Histo->Draw("HIST");
	//Dsmass_Histo->Draw("e1pSAME");	
	/*TLegend* leg_Dseta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Dseta->SetFillColor(kWhite);
	leg_Dseta->SetFillStyle(1001);
	leg_Dseta->SetBorderSize(0);
	//leg_Dsmass->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Dseta->AddEntry(Dseta_Histo,"Dseta","L");
	leg_Dseta->Draw();*/
	c2->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas2.png");
	//-------------------------------------------------------------------------
	TCanvas* c3 = new TCanvas("c3","Canvas 3 - behavior of the D",1200,600);
	c3->Divide(2);	
	c3->cd(1);
	//D0phi_Histo->SetLineColor(kBlue);
	D0phi_Histo->SetMarkerStyle(7);
	D0phi_Histo->Draw("HIST");
	//D0eta_Histo->Draw("e1pSAME");
	/*TLegend* leg_D0phi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_D0phi->SetFillColor(kWhite);
	leg_D0phi->SetFillStyle(1001);
	leg_D0phi->SetBorderSize(0);
	//leg_D0phi->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_D0phi->AddEntry(D0phi_Histo,"D0phi","L");
	leg_D0phi->Draw();*/
	//-------------------------------------------------------------------------
	c3->cd(2);	
	//Dsphi_Histo->SetLineColor(kBlue);
	Dsphi_Histo->SetMarkerStyle(7);
	Dsphi_Histo->Draw("HIST");
	//Dsphi_Histo->Draw("e1pSAME");	
	/*TLegend* leg_Dsphi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Dsphi->SetFillColor(kWhite);
	leg_Dsphi->SetFillStyle(1001);
	leg_Dsphi->SetBorderSize(0);
	//leg_Dsmass->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Dsphi->AddEntry(Dsphi_Histo,"Dsphi","L");
	leg_Dsphi->Draw();*/
	c3->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas3.png");

	if(Meson)
	{
	//=========================================================================	
	//Creating Canvas
	TCanvas* c4 = new TCanvas("c4","Canvas 4 - behavior of the SlowPion",1200,600);
	c4->Divide(2);
	c4->cd(1);
	//TrkSpt_Histo->SetLineColor(kBlue);
	TrkSpt_Histo->SetMarkerStyle(7);
	TrkSpt_Histo->Draw("HIST");
	//TrkSpt_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkSpt = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSpt->SetFillColor(kWhite);
	leg_TrkSpt->SetFillStyle(1001);
	leg_TrkSpt->SetBorderSize(0);
	//leg_TrkSpt->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSpt->AddEntry(TrkSpt_Histo,"TrkSpt","L");
	leg_TrkSpt->Draw();*/
	//-------------------------------------------------------------------------
	c4->cd(2);
	//TrkSnhits_Histo->SetLineColor(kBlue);
	TrkSnhits_Histo->SetMarkerStyle(7);
	TrkSnhits_Histo->Draw("HIST");
	//Dsphi_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkSnhits = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSnhits->SetFillColor(kWhite);
	leg_TrkSnhits->SetFillStyle(1001);
	leg_TrkSnhits->SetBorderSize(0);
	//leg_TrkSnhits->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSnhits->AddEntry(TrkSnhits_Histo,"TrkSnhits","L");
	leg_TrkSnhits->Draw();*/
	c4->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas4.png");
	//-------------------------------------------------------------------------
	TCanvas* c5 = new TCanvas("c5","Canvas 5 - behavior of the D",1200,600);
	c5->Divide(2);
	c5->cd(1);
	//TrkSdxy_Histo->SetLineColor(kBlue);
	TrkSdxy_Histo->SetMarkerStyle(7);
	TrkSdxy_Histo->Draw("HIST");
	//TrkSdxy_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkSdxy = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSdxy->SetFillColor(kWhite);
	leg_TrkSdxy->SetFillStyle(1001);
	leg_TrkSdxy->SetBorderSize(0);
	//leg_TrkSdxy->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSdxy->AddEntry(TrkSdxy_Histo,"TrkSdxy","L");
	leg_TrkSdxy->Draw();*/
	//-------------------------------------------------------------------------
	c5->cd(2);
	//TrkSdz_Histo->SetLineColor(kBlue);
	TrkSdz_Histo->SetMarkerStyle(7);
	TrkSdz_Histo->Draw("HIST");
	//TrkSdz_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkSdz = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSdz->SetFillColor(kWhite);
	leg_TrkSdz->SetFillStyle(1001);
	leg_TrkSdz->SetBorderSize(0);
	//leg_TrkSdz->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSdz->AddEntry(TrkSdz_Histo,"TrkSdz","L");
	leg_TrkSdz->Draw();*/
	c5->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas5.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c6 = new TCanvas("c6","Canvas 6 - behavior of the SlowPion",1200,600);
	c6->Divide(2);
	c6->cd(1);
	//TrkSeta_Histo->SetLineColor(kBlue);
	TrkSeta_Histo->SetMarkerStyle(7);
	TrkSeta_Histo->Draw("HIST");
	//TrkSeta_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkSeta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSeta->SetFillColor(kWhite);
	leg_TrkSeta->SetFillStyle(1001);
	leg_TrkSeta->SetBorderSize(0);
	//leg_TrkSeta->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSeta->AddEntry(TrkSeta_Histo,"TrkSeta","L");
	leg_TrkSeta->Draw();*/
	//-------------------------------------------------------------------------
	c6->cd(2);
	//TrkSphi_Histo->SetLineColor(kBlue);
	TrkSphi_Histo->SetMarkerStyle(7);
	//TrkSphi_Histo->SetMarkerStyle(21);
	//TrkSphi_Histo->SetStats(0);
	//TrkSphi_Histo->SetMarkerColor(kBlue);
	//TrkSphi_Histo->SetFillColor(kBlue);
	TrkSphi_Histo->Draw("HIST");
	//TrkSphi_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkSphi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSphi->SetFillColor(kWhite);
	leg_TrkSphi->SetFillStyle(1001);
	leg_TrkSphi->SetBorderSize(0);
	//leg_TrkSphi->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSphi->AddEntry(TrkSphi_Histo,"TrkSphi","L");
	leg_TrkSphi->Draw();*/
	c6->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas6.png");
	//-------------------------------------------------------------------------
	TCanvas* c7 = new TCanvas("c7","Canvas 7 - behavior of the D",1200,600);
	c7->Divide(2);
	c7->cd(1);
	//TrkSchi2_Histo->SetLineColor(kBlue);
	TrkSchi2_Histo->SetMarkerStyle(7);
	TrkSchi2_Histo->Draw("HIST");
	//TrkSchi2_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkSchi2 = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkSchi2->SetFillColor(kWhite);
	leg_TrkSchi2->SetFillStyle(1001);
	leg_TrkSchi2->SetBorderSize(0);
	//leg_TrkSchi2->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkSchi2->AddEntry(TrkSchi2_Histo,"TrkSchi2","L");
	leg_TrkSchi2->Draw();*/
	c7->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas7.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c8 = new TCanvas("c8","Canvas 8 - behavior of the Pion",1200,600);
	c8->Divide(2);
	c8->cd(1);
	//Trkpipt_Histo->SetLineColor(kBlue);
	Trkpipt_Histo->SetMarkerStyle(7);
	Trkpipt_Histo->Draw("HIST");
	//TrkSpt_Histo->Draw("e1pSAME");
	/*TLegend* leg_Trkpipt = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpipt->SetFillColor(kWhite);
	leg_Trkpipt->SetFillStyle(1001);
	leg_Trkpipt->SetBorderSize(0);
	//leg_Trkpipt->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpipt->AddEntry(Trkpipt_Histo,"Trkpipt","L");
	leg_Trkpipt->Draw();*/
	//-------------------------------------------------------------------------
	c8->cd(2);
	//Trkpinhits_Histo->SetLineColor(kBlue);
	Trkpinhits_Histo->SetMarkerStyle(7);
	Trkpinhits_Histo->Draw("HIST");
	//Trkpinhits_Histo->Draw("e1pSAME");
	/*TLegend* leg_Trkpinhits = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpinhits->SetFillColor(kWhite);
	leg_Trkpinhits->SetFillStyle(1001);
	leg_Trkpinhits->SetBorderSize(0);
	//leg_TrkSnhits->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpinhits->AddEntry(Trkpinhits_Histo,"Trkpinhits","L");
	leg_Trkpinhits->Draw();*/
	c8->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas8.png");
	//-------------------------------------------------------------------------
	//Creating Canvas
	TCanvas* c9 = new TCanvas("c9","Canvas 9 - behavior of the Pion",1200,600);
	c9->Divide(2);
	c9->cd(1);
	//Trkpidxy_Histo->SetLineColor(kBlue);
	Trkpidxy_Histo->SetMarkerStyle(7);
	Trkpidxy_Histo->Draw("HIST");
	//Trkpidxy_Histo->Draw("e1pSAME");
	/*TLegend* leg_Trkpidxy = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpidxy->SetFillColor(kWhite);
	leg_Trkpidxy->SetFillStyle(1001);
	leg_Trkpidxy->SetBorderSize(0);
	//leg_Trkpidxy->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpidxy->AddEntry(Trkpidxy_Histo,"Trkpidxy","L");
	leg_Trkpidxy->Draw();*/
	//-------------------------------------------------------------------------
	c9->cd(2);
	//Trkpidz_Histo->SetLineColor(kBlue);
	Trkpidz_Histo->SetMarkerStyle(7);
	Trkpidz_Histo->Draw("HIST");
	//Trkpidz_Histo->Draw("e1pSAME");
	/*TLegend* leg_Trkpidz = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpidz->SetFillColor(kWhite);
	leg_Trkpidz->SetFillStyle(1001);
	leg_Trkpidz->SetBorderSize(0);
	//leg_Trkpidz->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpidz->AddEntry(Trkpidz_Histo,"Trkpidz","L");
	leg_Trkpidz->Draw();*/
	c9->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas9.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c10 = new TCanvas("c10","Canvas 10 - behavior of the Pion",1200,600);
	c10->Divide(2);
	c10->cd(1);
	//Trkpieta_Histo->SetLineColor(kBlue);
	Trkpieta_Histo->SetMarkerStyle(7);
	Trkpieta_Histo->Draw("HIST");
	//Trkpieta_Histo->Draw("e1pSAME");
	/*TLegend* leg_Trkpieta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpieta->SetFillColor(kWhite);
	leg_Trkpieta->SetFillStyle(1001);
	leg_Trkpieta->SetBorderSize(0);
	//leg_Trkpieta->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpieta->AddEntry(Trkpieta_Histo,"Trkpieta","L");
	leg_Trkpieta->Draw();*/
	//-------------------------------------------------------------------------
	c10->cd(2);
	//Trkpiphi_Histo->SetLineColor(kBlue);
	Trkpiphi_Histo->SetMarkerStyle(7);
	Trkpiphi_Histo->Draw("HIST");
	//Trkpiphi_Histo->Draw("e1pSAME");
	/*TLegend* leg_Trkpiphi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpiphi->SetFillColor(kWhite);
	leg_Trkpiphi->SetFillStyle(1001);
	leg_Trkpiphi->SetBorderSize(0);
	//leg_Trkpiphi->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpiphi->AddEntry(Trkpiphi_Histo,"Trkpiphi","L");
	leg_Trkpiphi->Draw();*/
	c10->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas10.png");
	//-------------------------------------------------------------------------
	//Creating Canvas
	TCanvas* c11 = new TCanvas("c11","Canvas 11 - behavior of the Pion",1200,600);
	c11->Divide(2);
	c11->cd(1);
	//Trkpichi2_Histo->SetLineColor(kBlue);
	Trkpichi2_Histo->SetMarkerStyle(7);
	Trkpichi2_Histo->Draw("HIST");
	//Trkpichi2_Histo->Draw("e1pSAME");
	/*TLegend* leg_Trkpichi2 = new TLegend(0.82,0.5,0.95,0.65);
   	leg_Trkpichi2->SetFillColor(kWhite);
	leg_Trkpichi2->SetFillStyle(1001);
	leg_Trkpichi2->SetBorderSize(0);
	//leg_Trkpichi2->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_Trkpichi2->AddEntry(Trkpichi2_Histo,"Trkpichi2","L");
	leg_Trkpichi2->Draw();*/
	c11->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas11.png");

	//=========================================================================	
	//Creating Canvas
	TCanvas* c12 = new TCanvas("c12","Canvas 12 - behavior of the Kaon",1200,600);
	c12->Divide(2);
	c12->cd(1);
	//TrkKpt_Histo->SetLineColor(kBlue);
	TrkKpt_Histo->SetMarkerStyle(7);
	TrkKpt_Histo->Draw("HIST");
	//TrkSpt_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkKpt = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKpt->SetFillColor(kWhite);
	leg_TrkKpt->SetFillStyle(1001);
	leg_TrkKpt->SetBorderSize(0);
	//leg_TrkKpt->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKpt->AddEntry(TrkKpt_Histo,"TrkKpt","L");
	leg_TrkKpt->Draw();*/
	//-------------------------------------------------------------------------
	c12->cd(2);
	//TrkKnhits_Histo->SetLineColor(kBlue);
	TrkKnhits_Histo->SetMarkerStyle(7);
	TrkKnhits_Histo->Draw("HIST");
	//TrkKnhits_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkKnhits = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKnhits->SetFillColor(kWhite);
	leg_TrkKnhits->SetFillStyle(1001);
	leg_TrkKnhits->SetBorderSize(0);
	//leg_TrkSnhits->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKnhits->AddEntry(TrkKnhits_Histo,"TrkKnhits","L");
	leg_TrkKnhits->Draw();*/
	c12->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas12.png");
	//-------------------------------------------------------------------------
	//Creating Canvas
	TCanvas* c13 = new TCanvas("c13","Canvas 13 - behavior of the Pion",1200,600);
	c13->Divide(2);
	c13->cd(1);
	//TrkKdxy_Histo->SetLineColor(kBlue);
	TrkKdxy_Histo->SetMarkerStyle(7);
	TrkKdxy_Histo->Draw("HIST");
	//TrkKdxy_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkKdxy = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKdxy->SetFillColor(kWhite);
	leg_TrkKdxy->SetFillStyle(1001);
	leg_TrkKdxy->SetBorderSize(0);
	//leg_TrkKdxy->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKdxy->AddEntry(TrkKdxy_Histo,"TrkKdxy","L");
	leg_TrkKdxy->Draw();*/
	//-------------------------------------------------------------------------
	c13->cd(2);
	//TrkKdz_Histo->SetLineColor(kBlue);
	TrkKdz_Histo->SetMarkerStyle(7);
	TrkKdz_Histo->Draw("HIST");
	//TrkKdz_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkKdz = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKdz->SetFillColor(kWhite);
	leg_TrkKdz->SetFillStyle(1001);
	leg_TrkKdz->SetBorderSize(0);
	//leg_TrkKdz->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKdz->AddEntry(TrkKdz_Histo,"TrkKdz","L");
	leg_TrkKdz->Draw();*/
	c13->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas13.png");
	//=========================================================================	
	//Creating Canvas
	TCanvas* c14 = new TCanvas("c14","Canvas 14 - behavior of the Kaon",1200,600);
	c14->Divide(2);
	c14->cd(1);
	//TrkKeta_Histo->SetLineColor(kBlue);
	TrkKeta_Histo->SetMarkerStyle(7);
	TrkKeta_Histo->Draw("HIST");
	//TrkKeta_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkKeta = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKeta->SetFillColor(kWhite);
	leg_TrkKeta->SetFillStyle(1001);
	leg_TrkKeta->SetBorderSize(0);
	//leg_TrkKeta->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKeta->AddEntry(TrkKeta_Histo,"TrkKeta","L");
	leg_TrkKeta->Draw();*/
	//-------------------------------------------------------------------------
	c14->cd(2);
	//TrkKphi_Histo->SetLineColor(kBlue);
	TrkKphi_Histo->SetMarkerStyle(7);
	TrkKphi_Histo->Draw("HIST");
	//TrkKphi_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkKphi = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKphi->SetFillColor(kWhite);
	leg_TrkKphi->SetFillStyle(1001);
	leg_TrkKphi->SetBorderSize(0);
	//leg_TrkKphi->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKphi->AddEntry(TrkKphi_Histo,"TrkKphi","L");
	leg_TrkKphi->Draw();*/
	c14->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas14.png");
	//-------------------------------------------------------------------------
	//Creating Canvas
	TCanvas* c15 = new TCanvas("c15","Canvas 15 - behavior of the Pion",1200,600);
	c15->Divide(2);
	c15->cd(1);
	//TrkKchi2_Histo->SetLineColor(kBlue);
	TrkKchi2_Histo->SetMarkerStyle(7);
	TrkKchi2_Histo->Draw("HIST");
	//TrkKchi2_Histo->Draw("e1pSAME");
	/*TLegend* leg_TrkKchi2 = new TLegend(0.82,0.5,0.95,0.65);
   	leg_TrkKchi2->SetFillColor(kWhite);
	leg_TrkKchi2->SetFillStyle(1001);
	leg_TrkKchi2->SetBorderSize(0);
	//leg_TrkKchi2->AddEntry(h_leadingMuon_Pt2,"Data 2011","e1pSAME");
	leg_TrkKchi2->AddEntry(TrkKchi2_Histo,"TrkKchi2","L");
	leg_TrkKchi2->Draw();*/
	c15->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas15.png");
	}//End Meson Flag
	}//End VectorsHistograms Flag
}//end program

