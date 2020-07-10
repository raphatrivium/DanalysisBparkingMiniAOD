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
#include "TChain.h"
#ifndef ROOT_TLatex
#define ROOT_TLatex
#endif
#ifndef ROOTiosfwd
//#include "Riosfwd.h"
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
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "THStack.h"

//gSystem->AddIncludePath("/home/raphael/Desktop/analysis_2018");

using namespace std;

TChain *callTchain( const char* NAMETREE, const char* NAMEFILE);
TH1* GetTH1(const char* FILE, const char* hNAME);
TH1* makeTH1(const char* NAMETREE, const char* NAMEFILE , string BRANCH, string NBIN , string min , string max, string NAME);
TH1* makeTH1MC(const char* NAMETREE, const char* NAMEFILE, string BRANCH, string NBIN , string min , string max, string NAME, Color_t COLOR);
void DataMC(const char* NAMETREE_DATA, const char* NAMEFILE_DATA, const char* NAMETREE_MC, const char* NAMEFILE_MC, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR);
void PUDataMC(const char* NAMETREE_DATA, const char* NAMEFILE_DATA, const char* NAMEFILE_MC, string BRANCH, string NBIN, string min , string max, string NAME, string path, Color_t COLOR);
void PUDataMC2V(const char* NAMETREE_DATA, const char* NAMEFILE_DATA, const char* NAMEFILE_MC, const char* NAMEFILE_MC2,string BRANCH, string NBIN, string min , string max, string NAME, string path, Color_t COLOR);
void DataPUMCSigBkg(const char* NAMETREE_DATA, const char* NAMEFILE_DATA, const char* NAMEFILE_MC, const char* NAMEFILE_MC2,string BRANCH, string NBIN, string min , string max, string NAME, string path, Color_t COLOR);
void Drawhisto(const char* NAMETREE_DATA, const char* NAMEFILE_MC, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR);
void DrawhistoABS(const char* NAMETREE_DATA, const char* NAMEFILE_MC, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR);
void DrawhistoPow(const char* NAMETREE_DATA, const char* NAMEFILE_MC, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR);
void DrawhistoMClog(const char* NAMETREE, const char* NAMEFILE, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR);
void Drawhisto2D(string BRANCH1, string BRANCH2, string NBIN1 , string min1 , string max1, string NBIN2 , string min2 , string max2, string NAME, string path);
void TGraphAsymmErrorsRec(const char* FILE, const char* hNAME1, const char* hNAME2, const char* hTITLE , string path);


void makeHistograms(){
	gROOT->Reset();
	//call a file for a histogram style (Optional)
	//gROOT->LoadMacro("styleTDR.C");
	//gROOT->ProcessLine(".L styleTDR.C");  
	//styleTDR();

	//----------------------
	//FILE NAMES AND PATHS
	//----------------------
	//string path = "/home/raphael/cernbox/analysisB2019/canvas/"; //local to save the png histograms
	//string NameFileData = "/home/raphael/cernbox/testes/DMesonsHistograms_Data.root";
	//string NameFileMC = "/home/raphael/cernbox/testes/MC_DStarToD0Pi_D0KPi.root";

	string path = "canvas/DataMC/"; //local to save the histograms
	string path2 = "canvas/MC_DStarToD0Pi_D0KPi/";

	string NameFileData = "rootfit/BparkingDataset_Data.root";
	string NameFileMC = "rootfit/MC_DStarToD0Pi_D0KPi.root";
	string NameFilePUMC = "PUMC_DStarToD0Pi_D0KPi.root";
	string NameFileMinBiasPU = "PUMC_MinBias.root";

	//string path = "/eos/user/r/ragomesd/analysisB2019/canvas/"; //local to save the histograms
	//string path2 = "/eos/user/r/ragomesd/analysisB2019/canvas/MC_DStarToD0Pi_D0KPi/";
	//string NameFileData = "/eos/user/r/ragomesd/crab/haddteste/BparkingDataset_Data.root";
	//string NameFileMC = "/eos/user/r/ragomesd/crab/haddteste/MC_DStarToD0Pi_D0KPi_withoutTrigger.root";
	//string NameFileMC = "/eos/user/r/ragomesd/crab/haddteste/MC_DStarToD0Pi_D0KPi.root";
	//string NameFilePUMC = "/eos/user/r/ragomesd/crab/haddteste/PUMC_DStarToD0Pi_D0KPi.root";
	//string NameFileMinBiasPU = "/eos/user/r/ragomesd/crab/haddteste/PUMC_MinBias.root";
	//----------------------
	//TREE NAMES
	//----------------------
	string NameTreeDsRec = "t_analysis";  
	string NameTreeDsMC =  "t_DsMC";
	string NameTreeDsMat = "t_DsMatching";

	string NameTreeD0Rec = "t_D0analysis"; 
	string NameTreeD0MC =  "t_D0MC";
	string NameTreeD0Mat = "t_D0Matching";
	//------------------------------------
	//FLAGS
	//------------------------------------
	bool RecEff = true;
	bool flagMatching = true; //true - false
	bool mesonDstar = true;
	bool mesonD0 = false;
	
	TChain *chain = callTchain(NameTreeDsRec.c_str(),NameFileData.c_str());
	Long64_t nentries1 = chain->GetEntries(); //Reading Number of tree entries
	cout<< "Number of DATA entries : "<< nentries1 <<std::endl;
	TChain *chain2 = callTchain(NameTreeDsRec.c_str(),NameFileMC.c_str());
	Long64_t nentries2 = chain2->GetEntries(); //Reading Number of tree entries
	cout<< "Number of MC entries: "<< nentries2 <<std::endl;

	/*TChain *chain3 = callTchain(NameTreeDsMC.c_str(),NameFileData.c_str());
	TChain *chain4 = callTchain(NameTreeDsMC.c_str(),NameFileMC.c_str());

	TChain *chain5 = callTchain(NameTreeD0Rec.c_str(),NameFileData.c_str());
	TChain *chain6 = callTchain(NameTreeD0Rec.c_str(),NameFileMC.c_str());

	TChain *chain7 = callTchain(NameTreeD0MC.c_str(),NameFileData.c_str());
	TChain *chain8 = callTchain(NameTreeD0MC.c_str(),NameFileMC.c_str());*/


	//Long64_t nentries2 = chain2->GetEntries(); //Reading Number of tree entries
	//cout<< "Number of MC entries: "<< nentries2 <<std::endl;

	if (RecEff){
	TGraphAsymmErrorsRec(NameFileMC.c_str(),"hRecPt_Rec1", "hRecPt_MC", "Reconstruction Efficiency D*; p_{T} [GeV]; Efficiency", path2.c_str());
	TGraphAsymmErrorsRec(NameFileMC.c_str(),"hRecEta_Rec1", "hRecEta_MC", "Reconstruction Efficiency D*; #eta; Efficiency", path2.c_str());

	TGraphAsymmErrorsRec(NameFileMC.c_str(),"hRecD0Pt_Rec1", "hRecD0Pt_MC", "Reconstruction Efficiency D0; p_{T} [GeV]; Efficiency", path2.c_str());
	TGraphAsymmErrorsRec(NameFileMC.c_str(),"hRecD0Eta_Rec1", "hRecD0Eta_MC", "Reconstruction Efficiency D0; #eta; Efficiency", path2.c_str());
	}
	//---------------------------------------------------------------------------------------------
	string Branch = ""; string nbin = ""; string nmin = ""; string nmax = ""; string hTitle = "";
	//--------------------------------------------------------------------------------
	if(flagMatching)
	{
	Branch = "deltaRDsMatching"; nbin = "100"; nmin = "0"; nmax = "0.06"; hTitle = "#DeltaR Ds ; #DeltaR ; Events ";
	DrawhistoMClog(NameTreeDsMat.c_str(), NameFileMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path2.c_str(), 0);
	Branch = "deltaRD0Matching"; nbin = "100";  nmin = "0";  nmax = "0.06"; hTitle = "#DeltaR D0 ; #DeltaR ; Events ";
	DrawhistoMClog(NameTreeD0Mat.c_str(), NameFileMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path2.c_str(), 0);
	}
	//--------------------------------------------------------------------------------
	if(mesonDstar){
	Branch = "Dsmass";  nbin = "100";  nmin = "1.93";  nmax = "2.1";	hTitle = "Invariant Mass of the From D* ; Mass [GeV] ; Events ";
	DataMC(NameTreeDsRec.c_str(), NameFileData.c_str(), NameTreeDsRec.c_str(), NameFileMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	PUDataMC(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	PUDataMC2V(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	
	Branch = "Dspt";  nbin = "100";  nmin = "0";  nmax = "30"; hTitle = "p_{T} distribuition of the D* ; p_{T} [GeV] ; Events ";
	DataMC(NameTreeDsRec.c_str(), NameFileData.c_str(), NameTreeDsRec.c_str(), NameFileMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	PUDataMC(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	PUDataMC2V(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);

	Branch = "D0Kpipt";  nbin = "100";  nmin = "0";  nmax = "30"; hTitle = "p_{T} distribuition of the D0 ; p_{T} [GeV] ; Events ";
	DataPUMCSigBkg(NameTreeD0Rec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	/*Branch = "Dseta";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "Pseudo-rapidity distribuition of the D* ; #eta ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "Dsphi";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "#Phi distribuition of the D* ; #Phi ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "D0mass";  nbin = "100";  nmin = "1.77";  nmax = "1.95";	hTitle = "Invariant Mass of the D0(From D*) ; Mass [GeV] ; Events";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "D0pt";  nbin = "100";  nmin = "0";  nmax = "30"; hTitle = "pT distribuition of the D0(From D*) ; p_{T} [GeV] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	
	Branch = "D0eta";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "Pseudo-rapidity distribuition of the D0(From D*) ; #eta ; Events";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "D0phi";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "#Phi distribuition of the D0(From D*) ; #Phi ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkSpt";  nbin = "100";  nmin = "0";  nmax = "10"; hTitle = "pt distribuition of the Kaon ; p_{T} [GeV] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkSnhits";  nbin = "100";  nmin = "0";  nmax = "50"; hTitle = "nhits distribuition of the SlowPion ; Number of hits ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkSdxy";  nbin = "100";  nmin = "-0.1";  nmax = "0.1"; hTitle = "dxy distribuition of the SlowPion ; dx [cm] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkSdz";  nbin = "100";  nmin = "-0.1";  nmax = "0.1"; hTitle = "dz distribuition of the SlowPion ; dz [cm] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkSeta";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "eta distribuition of the SlowPion ; #eta ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkSphi";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "Phi distribuition of the SlowPion ; #Phi ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkSchi2";  nbin = "100";  nmin = "0";  nmax = "4"; hTitle = "Chi2 distribuition of the SlowPion ; #Chi^{2} ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "Trkpipt";  nbin = "100";  nmin = "0";  nmax = "30"; hTitle = "pt distribuition of the Pion ; p_{T} [GeV] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "Trkpinhits";  nbin = "100";  nmin = "0";  nmax = "50"; hTitle = "nhits distribuition of the Pion ; Number of hits ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "Trkpidxy";  nbin = "100";  nmin = "-0.1";  nmax = "0.1"; hTitle = "dxy distribuition of the Pion ; dx [cm] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "Trkpidz";  nbin = "100";  nmin = "-0.1";  nmax = "0.1"; hTitle = "dz distribuition of the Pion ; dz [cm] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "Trkpieta";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "eta distribuition of the Pion ; #eta ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "Trkpiphi";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "phi distribuition of the Pion ; #Phi ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "Trkpichi2";  nbin = "100";  nmin = "0";  nmax = "4"; hTitle = "Chi2 distribuition of the Pion ; #Chi^{2} ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkKpt";  nbin = "100";  nmin = "0";  nmax = "30"; hTitle = "pt distribuition of the Kaon ; p_{T} [GeV] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkKnhits";  nbin = "100";  nmin = "0";  nmax = "50"; hTitle = "Invariant Mass of the D0(From D*) ; Mass [GeV] ; Events";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkKdxy";  nbin = "100";  nmin = "-0.1";  nmax = "0.1"; hTitle = "dxy distribuition of the Kaon ; dx [cm] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkKdz";  nbin = "100";  nmin = "-0.1";  nmax = "0.1"; hTitle = "dz distribuition of the Kaon ; dz [cm] ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkKeta";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "eta distribuition of the Kaon ; #eta ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkKphi";  nbin = "100";  nmin = "-4";  nmax = "4"; hTitle = "phi distribuition of the Kaon ; #Phi ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	Branch = "TrkKchi2";  nbin = "100";  nmin = "0";  nmax = "4"; hTitle = "Chi2 distribuition of the Kaon ; #Chi^{2} ; Events ";
	DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);*/
	//Branch = "D0fromDSs3D";  nbin = "100";  nmin = "0";  nmax = "5"; hTitle = "3D Significance  of D0(From D*) ; Significance ; Events ";
	//DataPUMCSigBkg(NameTreeDsRec.c_str(), NameFileData.c_str(), NameFileMinBiasPU.c_str(), NameFilePUMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kRed);
	//Branch = "Anglephi";  nbin = "100";  nmin = "0.985";  nmax = "1.01"; hTitle = "#phi pointing of D0(From D*) ; #phi pointing ; Events ";
	//DataMC(NameTreeDsRec.c_str(), NameFileData.c_str(), NameTreeDsRec.c_str(), NameFileMC.c_str(), Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(), hTitle.c_str(), path.c_str(), kGreen);
	//Branch = "D0_VtxProb";  nbin = "100";  nmin = "0";  nmax = "1";
	//Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"D0_VtxProb ; Probability ; Events ", path.c_str(),kRed);
	}
	//--------------------------------------------------------------------------------
	/*if(mesonD0){
	//Variables D0 Prompt
	Branch = "D0Kpipt";  nbin = "100";  nmin = "0";  nmax = "30";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"p_{T} distribuition of the D0Prompt ; #eta ; Events ", path.c_str(),kRed);
	Branch = "D0Kpieta";  nbin = "100";  nmin = "-4";  nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"eta distribuition of the D0Prompt ; #eta ; Events ", path.c_str(),kRed);
	Branch = "D0Kpiphi";  nbin = "100";  nmin = "-4";  nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"phi distribuition of the D0Prompt ; #Phi ; Events ", path.c_str(),kRed);
	Branch = "TrkD0Kdxy";  nbin = "100";  nmin = "-0.1";  nmax = "0.1";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dx distribuition of the Kaon(D0Prompt) ; dz [cm] ; Events  ", path.c_str(),kRed);
	Branch = "TrkD0pidxy";  nbin = "100";  nmin = "-0.1";  nmax = "0.1";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dx distribuition of the pion(D0Prompt) ; dz [cm] ; Events  ", path.c_str(),kRed);
	Branch = "TrkD0Kdz";  nbin = "100";  nmin = "-0.1";  nmax = "0.1";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dz distribuition of the Kaon(D0Prompt) ; dz [cm] ; Events ", path.c_str(),kRed);
	Branch = "TrkD0pidz";  nbin = "100";  nmin = "-0.1";  nmax = "0.1";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dz distribuition of the pion(D0Prompt) ; dz [cm] ; Events ", path.c_str(),kRed);
	Branch = "TrkD0Kchi2";  nbin = "100";  nmin = "0";  nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"Chi2 distribuition of the Kaon(D0Prompt) ; #Chi^{2} ; Events ", path.c_str(),kRed);
	Branch = "TrkD0pichi2";  nbin = "100";  nmin = "0";  nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"Chi2 distribuition of the pion(D0Prompt) ; #Chi^{2} ; Events ", path.c_str(),kRed);
	Branch = "TrkD0Kpt";  nbin = "100";  nmin = "0";  nmax = "30";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"p_{T} distribuition of the Kaon(D0Prompt) ; pT [GeV] ; Events ", path.c_str(),kRed);
	Branch = "TrkD0pipt";  nbin = "100";  nmin = "0";  nmax = "30";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"p_{T} distribuition of the pion(D0Prompt) ; pT [GeV] ; Events ", path.c_str(),kRed);
	Branch = "TrkD0Keta";  nbin = "100";  nmin = "-4";  nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"eta distribuition of the Kaon(D0Prompt) ; #eta ; Events  ", path.c_str(),kRed);
	Branch = "TrkD0pieta";  nbin = "100";  nmin = "-4";  nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"eta distribuition of the pion(D0Prompt) ; #eta ; Events  ", path.c_str(),kRed);
	Branch = "TrkD0kphi"; nbin = "100";  nmin = "-4";  nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"phi distribuition of the kaon(D0Prompt) ; #Phi ; Events ", path.c_str(),kRed);
	Branch = "TrkD0piphi";  nbin = "100";  nmin = "-4";  nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"phi distribuition of the pion(D0Prompt) ; #Phi ; Events", path.c_str(),kRed);
	Branch = "TrkD0Knhits";  nbin = "100";  nmin = "0";  nmax = "50";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"nhits distribuition of the Kaon(D0Prompt) ; Number of hits ; Events  ", path.c_str(),kRed);
	Branch = "TrkD0pinhits";  nbin = "100";  nmin = "0";  nmax = "50";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"nhits distribuition of the pion(D0Prompt) ; Number of hits ; Events  ", path.c_str(),kRed);
	}*/

	cout << "----------------------" << endl;
	cout << "E N D   P R O G R A M" << endl;
	cout << "----------------------" << endl;
	
}//end program
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TChain *callTchain( const char* NAMETREE, const char* NAMEFILE){

	TChain* chain = new TChain(NAMETREE);
	chain->Add(NAMEFILE);
	return chain;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH1* GetTH1(const char* FILE, const char* hNAME)
{	//Get the histogram from a tree
	TFile *file = new TFile(FILE);
	TH1D *h1 = (TH1D*)file->Get(hNAME);
	return h1;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH1* makeTH1(const char* NAMETREE, const char* NAMEFILE , string BRANCH, string NBIN , string min , string max, string NAME)
{
	//To make histograms from the data file
	//TChain chain("t_analysis"); chain.Add("DMesonsHistograms_Data.root");
	TChain *chain1 = callTchain(NAMETREE,NAMEFILE);
	chain1->Draw(Form("%s>>HistoDATA(%s,%s,%s)",BRANCH.c_str(), NBIN.c_str(), min.c_str(), max.c_str()));
	TH1 *hhData = (TH1*)gDirectory->Get("HistoDATA");
	hhData->SetTitle(NAME.c_str());
	hhData->SetMarkerStyle(21);
	hhData->Sumw2();
	//hh->Draw("HIST");
   return hhData ;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH1* makeTH1MC(const char* NAMETREE, const char* NAMEFILE, string BRANCH, string NBIN , string min , string max, string NAME, Color_t COLOR)
{
	//To make histograms from the MC file
	TChain *chain2 = callTchain(NAMETREE, NAMEFILE);
	chain2->Draw(Form("%s>>HistoMC(%s,%s,%s)",BRANCH.c_str(), NBIN.c_str(), min.c_str(), max.c_str()));
	TH1 *hhMC = (TH1*)gDirectory->Get("HistoMC");
	hhMC->SetTitle(NAME.c_str());
	hhMC->SetFillStyle(4050); hhMC->SetFillColor(46);
	hhMC->SetLineColor(46);   hhMC->SetLineWidth(1);
	hhMC->SetYTitle("A.U");   hhMC->SetFillColorAlpha(COLOR,0.35);
	//hhMC->SetFillColor(COLOR);
	hhMC->SetFillStyle(4050);
	hhMC->Sumw2(); return hhMC ;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH1* PUmakeTH1MC(const char* NAMEFILE, string BRANCH, string NAME, Color_t COLOR)
{
	//To make histograms from the MC file with pileup weight
	TFile *file = new TFile(NAMEFILE);
	TH1D *hhMC = (TH1D*)file->Get(Form("h%s", BRANCH.c_str()));
	hhMC->SetTitle(NAME.c_str());
	hhMC->SetFillStyle(4050); hhMC->SetFillColor(46);
	hhMC->SetLineColor(46);   hhMC->SetLineWidth(1);
	hhMC->SetYTitle("A.U");   hhMC->SetFillColorAlpha(COLOR,0.35);
	hhMC->SetFillStyle(4050);
	hhMC->Sumw2(); return hhMC ;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void DataMC(const char* NAMETREE_DATA, const char* NAMEFILE_DATA, const char* NAMETREE_MC, const char* NAMEFILE_MC, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR)
{  //To put the Data and MC in the same histo with scale
	TCanvas* canvas = new TCanvas("c0","Canvas 0 - behavior of the D",1200,600);
	TH1 *Histodata = makeTH1(NAMETREE_DATA, NAMEFILE_DATA, BRANCH.c_str(), NBIN.c_str(), min.c_str(), max, NAME);
	TH1 *HistoMC = makeTH1MC(NAMETREE_MC, NAMEFILE_MC, BRANCH.c_str(), NBIN.c_str(), min.c_str(), max, NAME, COLOR);
	Histodata->SetStats(0); HistoMC->SetStats(0);
	//HistoMC->Scale(0.990);
	//HistoMC->Scale(HistoMC->Integral("width"));	
	//double scalevar=((Histodata->Integral())/(HistoMC->Integral()));
	//HistoMC->Scale(scalevar);
	//Double_t scale1 = 1./(Histodata->Integral()); Histodata->Scale(scale1);
	//Double_t scale2 = 1./(HistoMC->Integral()); HistoMC->Scale(scale2);
	HistoMC->Scale(Histodata->Integral()/HistoMC->Integral());
	HistoMC->Draw("Histo");	Histodata->Draw("e1pSAME");
	//D0massHistoMC->DrawNormalized();
	TLegend* leg = new TLegend(0.72,0.76,0.90,0.90);
   //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   leg->AddEntry(Histodata,"Data Bparking 2018","P");	
	leg->AddEntry(HistoMC,"MC","f");
   leg->Draw();

	canvas->SaveAs(Form("%s%s.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PUDataMC(const char* NAMETREE_DATA, const char* NAMEFILE_DATA, const char* NAMEFILE_MC, string BRANCH, string NBIN, string min , string max, string NAME, string path, Color_t COLOR)
{  //To put the Data and MC in the same histo with scale
	TCanvas* canvas = new TCanvas("c0","Canvas 0 - behavior of the D",1200,600);
	TH1 *Histodata = makeTH1(NAMETREE_DATA, NAMEFILE_DATA, BRANCH.c_str(), NBIN.c_str(), min.c_str(), max, NAME);
	TH1 *HistoMC = PUmakeTH1MC(NAMEFILE_MC, BRANCH.c_str(), NAME, COLOR); 
	Histodata->SetStats(0); HistoMC->SetStats(0);
	//double MCweight = ((5.112*pow(10,8)*0.000647218966)/10440479);
	double MCweight = ((7.842*pow(10,10)*0.000647218966)/9986000); //MinBias
	//HistoMC->Scale(MCweight);
	HistoMC->Scale(Histodata->Integral()/HistoMC->Integral());
	std::cout << "MCweight: " << MCweight << std::endl;
	//Histodata->Draw("e1p"); HistoMC->Draw("HistoSAME");
	HistoMC->Draw("Histo"); Histodata->Draw("e1pSAME");
	TLegend* leg = new TLegend(0.72,0.76,0.90,0.90);
   leg->AddEntry(Histodata,"Data Bparking 2018","P");	
	leg->AddEntry(HistoMC,"MC","f");
   leg->Draw();
	canvas->SaveAs(Form("%s%sPU.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PUDataMC2V(const char* NAMETREE_DATA, const char* NAMEFILE_DATA, const char* NAMEFILE_MC, const char* NAMEFILE_MC2,string BRANCH, string NBIN, string min , string max, string NAME, string path, Color_t COLOR)
{  //To put the Data and MC in the same histo with scale
	TCanvas* canvas = new TCanvas("c0","Canvas 0 - behavior of the D",1200,600);
	TH1 *Histodata = makeTH1(NAMETREE_DATA, NAMEFILE_DATA, BRANCH.c_str(), NBIN.c_str(), min.c_str(), max, NAME);
	TH1 *HistoMC = PUmakeTH1MC(NAMEFILE_MC, BRANCH.c_str(), NAME, COLOR); 
	TH1 *HistoMC2 = PUmakeTH1MC(NAMEFILE_MC2, BRANCH.c_str(), NAME, kGreen);
	Histodata->SetStats(0); HistoMC->SetStats(0); HistoMC2->SetStats(0);
	double MCweight = ((7.842*pow(10,10)*647.218966)/9986000); //MinBias
	//HistoMC->Scale(1./HistoMC->Integral()); HistoMC2->Scale(1./HistoMC2->Integral()); Histodata->Scale(1./Histodata->Integral());
	//HistoMC->Scale(Histodata->GetEntries()/HistoMC->GetEntries()); 
	//HistoMC2->Scale(Histodata->GetEntries()/HistoMC2->GetEntries()); //Histodata->Scale(Histodata->Integral());
	//std::cout << "MCweight: " << MCweight << std::endl;
	MCweight = ((5.112*pow(10,8)*647.218966)/10440479);
	HistoMC->Scale(Histodata->Integral()/HistoMC->Integral());	
	HistoMC2->Scale(Histodata->Integral()/HistoMC2->Integral());	
	//std::cout << "MCweight: " << MCweight << std::endl;
	//Histodata->Draw("e1p"); HistoMC->Draw("HistoSAME");
	HistoMC2->Draw("Histo"); HistoMC->Draw("HistoSAME"); Histodata->Draw("e1pSAME"); 
	//Histodata->Draw("e1p"); HistoMC->Draw("HistoSAME"); HistoMC2->Draw("HistoSAME");
	TLegend* leg = new TLegend(0.72,0.76,0.90,0.90);
   leg->AddEntry(Histodata,"Data Bparking 2018","P");	
	leg->AddEntry(HistoMC,"MC MinBias","f");
	leg->AddEntry(HistoMC2,"MC Signal","f");
   leg->Draw();
	canvas->SaveAs(Form("%s%sPU2.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void DataPUMCSigBkg(const char* NAMETREE_DATA, const char* NAMEFILE_DATA, const char* NAMEFILE_MC, const char* NAMEFILE_MC2, string BRANCH, string NBIN, string min , string max, string NAME, string path, Color_t COLOR)
{  //To put the Data and MC in the same histo with scale
	TCanvas* canvas = new TCanvas("c0","Canvas 0 - behavior of the D",1200,600);
	TH1 *Histodata = makeTH1(NAMETREE_DATA, NAMEFILE_DATA, BRANCH.c_str(), NBIN.c_str(), min.c_str(), max, NAME);
	TH1 *HistoMC = PUmakeTH1MC(NAMEFILE_MC, BRANCH.c_str(), NAME, COLOR); 
	TH1 *HistoMC2 = PUmakeTH1MC(NAMEFILE_MC2, BRANCH.c_str(), NAME, kGreen);
	Histodata->SetStats(0); HistoMC->SetStats(0); HistoMC2->SetStats(0);
	//double MCweight = ((7.842*pow(10,10)*647.218966)/9986000); //MinBias
	//HistoMC->Scale(20.);
	HistoMC->Scale(Histodata->Integral()/HistoMC->Integral());	
	//MCweight = ((5.112*pow(10,8)*647.218966)/10440479);
	//HistoMC2->Scale(0.15);
	HistoMC2->Scale(Histodata->Integral()/HistoMC2->Integral());	
	THStack hs("hs","hs");
	hs.SetTitle(NAME.c_str());
	hs.Add(HistoMC); hs.Add(HistoMC2);
	hs.Draw("Histo"); Histodata->Draw("e1pSAME");
	//Histodata->Draw("e1p"); hs.Draw("HistoSAME");
	TLegend* leg = new TLegend(0.72,0.76,0.90,0.90);
   leg->AddEntry(Histodata,"Data Bparking 2018","P");	
	leg->AddEntry(HistoMC,"MC MinBias","f");
	leg->AddEntry(HistoMC2,"MC Signal","f");
   leg->Draw();
	canvas->SaveAs(Form("%s%sSigBkg.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Drawhisto(const char* NAMETREE_DATA, const char* NAMEFILE_MC, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR)
{	//Only to draw histograms
	TChain *chain = callTchain(NAMETREE_DATA,NAMEFILE_MC);//"/home/raphael/cernbox/testes/DMesonsHistograms_Data.root"
	TCanvas* canvas = new TCanvas("canvas","",1200,600);
	chain->Draw(Form("%s>>Histo(%s,%s,%s)",BRANCH.c_str(), NBIN.c_str(), min.c_str(), max.c_str()));
	TH1 *Teste = (TH1*)gDirectory->Get("Histo");
	Teste->SetTitle(NAME.c_str());
	Teste->SetMarkerStyle(7);
	Teste->SetFillColor(COLOR);
	Teste->Sumw2();
	Teste->Draw("HIST");	
	canvas->SaveAs(Form("%s%s.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void DrawhistoABS(const char* NAMETREE_DATA, const char* NAMEFILE_MC, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR)
{	//Only to draw histograms
	TChain *chain = callTchain(NAMETREE_DATA,NAMEFILE_MC);//"/home/raphael/cernbox/testes/DMesonsHistograms_Data.root"
	TCanvas* canvas = new TCanvas("canvas","",1200,600);
	chain->Draw(Form("abs(%s)>>Histo(%s,%s,%s)",BRANCH.c_str(), NBIN.c_str(), min.c_str(), max.c_str()));
	TH1 *Teste = (TH1*)gDirectory->Get("Histo");
	Teste->SetTitle(NAME.c_str());
	Teste->SetMarkerStyle(7);
	Teste->SetFillColor(COLOR);
	Teste->Sumw2();
	Teste->Draw("HIST");	
	canvas->SaveAs(Form("%s%sABS.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void DrawhistoPow(const char* NAMETREE_DATA, const char* NAMEFILE_MC, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR)
{	//Only to draw histograms
	TChain *chain = callTchain(NAMETREE_DATA,NAMEFILE_MC);//"/home/raphael/cernbox/testes/DMesonsHistograms_Data.root"
	TCanvas* canvas = new TCanvas("canvas","",1200,600);
	chain->Draw(Form("pow(%s,2)>>Histo(%s,%s,%s)",BRANCH.c_str(), NBIN.c_str(), min.c_str(), max.c_str()));
	TH1 *Teste = (TH1*)gDirectory->Get("Histo");
	Teste->SetTitle(NAME.c_str());
	Teste->SetMarkerStyle(7);
	Teste->SetFillColor(COLOR);
	Teste->Sumw2();
	Teste->Draw("HIST");	
	canvas->SaveAs(Form("%s%sPow.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void DrawhistoMClog(const char* NAMETREE, const char* NAMEFILE, string BRANCH, string NBIN , string min , string max, string NAME, string path, Color_t COLOR)
{	//Only to draw histograms
	TChain *chain = callTchain(NAMETREE, NAMEFILE); //"/home/raphael/cernbox/testes/MC_DStarToD0Pi_D0KPi.root"
	TCanvas* canvas = new TCanvas("canvas","",1200,600);
	chain->Draw(Form("%s>>Histo(%s,%s,%s)",BRANCH.c_str(), NBIN.c_str(), min.c_str(), max.c_str()));
	TH1 *Teste = (TH1*)gDirectory->Get("Histo");
	Teste->SetTitle(NAME.c_str());
	Teste->SetStats(0); Teste->SetMarkerStyle(7);
	Teste->SetFillColor(COLOR); Teste->Sumw2();
	Teste->Draw("HIST"); gPad->SetLogy();
	canvas->SaveAs(Form("%s%s.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Drawhisto2D(string BRANCH1, string BRANCH2, string NBIN1 , string min1 , string max1, string NBIN2 , string min2 , string max2, string NAME, string path)
{	//To Draw 2D histograms using 
	TChain *chain = callTchain("t_analysis","/home/raphael/cernbox/testes/DMesonsHistograms_Data.root");
	TCanvas* canvas = new TCanvas("canvas","",1200,600);
	chain->Draw(Form("%s:%s>>Histo(%s,%s,%s,%s,%s,%s)",BRANCH1.c_str(), BRANCH2.c_str(), NBIN1.c_str(), min1.c_str(), max1.c_str(),
 																													 NBIN2.c_str(), min2.c_str(), max2.c_str()));
	TH2 *Teste = (TH2*)gDirectory->Get("Histo");
	Teste->SetTitle(NAME.c_str());
	Teste->Draw("COLZ");	//COL,COLZ,CONT0,CONTZ,CONT4COLZ
	canvas->SaveAs(Form("%sTH2%s%s.png", path.c_str(), BRANCH1.c_str(), BRANCH2.c_str()));
	delete canvas;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void TGraphAsymmErrorsRec(const char* FILE, const char* hNAME1, const char* hNAME2, const char* hTITLE, string path)
{	//Make a TGraphAsymmErrors usinh histograms as input
	TH1 * h1 = GetTH1(FILE, hNAME1);
	TH1 * h2 = GetTH1(FILE, hNAME2);
	if( TEfficiency::CheckConsistency(*h1,*h2) )
	{	TCanvas* canvas = new TCanvas("canvas","",1200,600);
		TGraphAsymmErrors *gr = new TGraphAsymmErrors(h1,h2);
		gr->SetTitle(hTITLE);
		gr->SetMarkerColor(4);
		gr->SetMarkerStyle(21);
		gr->Draw("AP");
		canvas->SaveAs(Form("%sTGraphAsymmErrors%s.pdf", path.c_str(), hNAME1));
		gr->GetErrorY(1);
		std::cout<< "gr->GetErrorY(0): " << gr->GetErrorY(0) << std::endl;
		std::cout<< "gr->GetErrorY(1): " << gr->GetErrorY(1) << std::endl;
		std::cout<< "gr->GetErrorY(2): " << gr->GetErrorY(2) << std::endl;
		std::cout<< "gr->GetErrorY(3): " << gr->GetErrorY(3) << std::endl;
		std::cout<< "gr->GetErrorY(4): " << gr->GetErrorY(4) << std::endl;
		std::cout<< "gr->GetErrorY(5): " << gr->GetErrorY(5) << std::endl;
		std::cout<< "gr->GetErrorY(6): " << gr->GetErrorY(6) << std::endl;
		std::cout<< "gr->GetErrorY(7): " << gr->GetErrorY(7) << std::endl;
		std::cout<< "gr->GetErrorY(8): " << gr->GetErrorY(8) << std::endl;
		std::cout<< "gr->GetErrorY(9): " << gr->GetErrorY(9) << std::endl;
		std::cout<< "gr->GetErrorY(10): " << gr->GetErrorY(10) << std::endl;
		std::cout<< "gr->GetErrorY(11): " << gr->GetErrorY(11) << std::endl;
		std::cout<< "gr->GetErrorY(12): " << gr->GetErrorY(12) << std::endl;
		double Eff = (gr->Integral(0,-1))/(gr->Integral());
		std::cout << "Eff: " << Eff << std::endl;
		delete canvas;
	}
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
