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
#include "TH1.h"
#include "TRandom.h"
#include <vector>
//#ifndef ROOT_TLatex
//#define ROOT_TLatex
//#ifndef ROOTiosfwd
#include "Riosfwd.h"
//#endif
//#ifndef ROOT_TText
#include "TText.h"
//#endif
//#ifndef ROOT_TAttLine
#include "TAttLine.h"
//#endif
//#include "RooCBShape.h"  //Crystal Ball
//#ifndef __CINT__
//#include "RooGlobalFunc.h"
//#endif
//#include "RooRealVar.h"
//#include "RooDataSet.h"
//#include "RooGaussian.h"
#include "TCanvas.h"
//#include "RooPlot.h"
#include "TAxis.h"
//using namespace RooFit ;
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TChain.h"
#include <time.h>
#include "TLegend.h"

//gSystem->AddIncludePath("/home/raphael/Desktop/analysis_2018");

using namespace std;
//bool IsOdd (string i) {return (i == 'HLT_Mu9_IP6');}

//bool IsOdd(string i) { return i = 2;}

//bool FindTrigger(string i) { return i = 2;}

//(s.find("Hello") == 0)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TChain *callOneFile(const char* TITLES){

	TChain* chain = new TChain("analysis/data");
	chain->Add(TITLES);
	return chain;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TChain *callTchain(const char* LOCAL){
	//READING ALL THE ROOT FILE IN THE LOCAL FOLDER USING TCHAIN
	const char *dirname = LOCAL; //Directory of the files 
	const char *ext=".root"; //extenson of files you want add	
	int added = 0;
	TChain* chain = new TChain("analysis/data"); 
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
	if (files) 
	{
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) 
		{
         fname = file->GetName();
			if (fname.EndsWith(ext))
			{	cout << "File: "<< fname << endl; 
				chain->Add(LOCAL+fname);
				added++;}	       
     }
   }
	cout << "======== " << added << " Files added"<< endl;
	return chain;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH1* makeTH1(const char* name, Double_t BIN, Double_t NMIN, Double_t NMAX, const char* TITLES, Color_t COLOR) 
{
   //Create a TH1D (one dimension) histogram
   TH1D* hh = new TH1D(name, name, BIN, NMIN, NMAX) ;
	hh->SetTitle(TITLES);
	//DsmassHisto->SetName("DsmassHisto");
	//hh->SetTitle(Form("%s", TITLES.c_str()));
	hh->SetMarkerStyle(7);
	hh->SetFillColor(COLOR); //kRed
	hh->Sumw2();
	//hh->Draw(TYPE); const char* TYPE
   return hh ;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//TBranch* makeBranch(TTree *myTree, const char* name, std::vector<double>* VECTOR) 
//{
//	TBranch *branch = myTree->GetBranch(name);
//	branch->SetAddress(&VECTOR);
//	return branch ;
//}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void analysisB2019_OfflineCuts()
{
	clock_t tStart = clock();

	string Dataset = "1"; 
	//Datasets: BparkingDataset1 - BparkingDataset2 - BparkingDataset3
	//Datasets: BparkingDataset4 - BparkingDataset5 - BparkingDataset6
	//Datasets: MC_DStarToD0Pi_D0KPi

	//Save the file
	string path = "/eos/user/r/ragomesd/crab/haddteste/";
	double SigCutD0fromDs = 1.0; //3.
	double SigCutD0 = 1.0; //3.

	bool debug = false;
	//eos/user/r/ragomesd/crab
		
	TChain* t1 = callTchain("/eos/user/r/ragomesd/crab/DStarToD0Pi_D0KPi_NoMuFilter_TuneCP5_13TeV-pythia8-evtgen/Ds_MC_Run2018A_Official_Tracks_HLT_Mu9_IP6/191213_204234/0000/");
	//TChain* t1 = callOneFile("/eos/user/r/ragomesd/testes/D0DstarDataBparking.root");
	//TChain* t1 = new TChain("analysis/data"); t1->Add("/eos/user/r/ragomesd/crab/ParkingBPH1/ParkingBPH_Run2018A_MINIAOD/191212_023520/0000/D0DstarDataBparking_35.root");
	
	//-------Reading the root file-------------------------------------	
	//TFile *f1 = new TFile("D0DstarDataBparking_7.root");
	//TTree *t1 = (TTree*)f1->Get("analysis/data");
	
	//-------Create root file-------------------------------------
	//TFile f_analysis("DMesonsHistograms.root","recreate"); //Creates root file
	TFile f_analysis(Form("%sDMesonsHistograms%s.root", path.c_str(), Dataset.c_str() ),"recreate"); //Creates root file
	TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	TTree t_DstarWrongCombination("t_DstarWrongCombination","t_DstarWrongCombination"); //Creates a Tree
	TTree t_D0analysis("t_D0analysis","t_D0analysis"); //Creates a Tree

	//It must be equal to "0" and not "0.", even it is a double. Possible segmental error if you compile it with c++.
	std::vector<string>* NameOfFiredTriggers = 0;

	//----------------------------------------
	//For D* 
	//-----------------------------------------
	std::vector<double>* D0mass = 0;
	std::vector<double>* D0eta = 0;
	std::vector<double>* D0phi = 0;
	std::vector<double>* D0pt = 0;

	std::vector<double>* Dsmass = 0;
	std::vector<double>* Dseta = 0;
	std::vector<double>* Dsphi = 0;
	std::vector<double>* Dspt = 0;

	//std::vector<double>* TrkScharge = 0;
	std::vector<double>* D0fromDSsXY = 0;
	std::vector<double>* D0fromDSs3D = 0;
	std::vector<double>* Anglephi = 0;
	std::vector<double>* D0_VtxProb = 0;

	std::vector<double>* TrkKpt = 0;
	std::vector<double>* TrkKnhits = 0;
	std::vector<double>* TrkKchi2 = 0;
	std::vector<double>* TrkKdxy = 0;
	std::vector<double>* TrkKdz = 0;
	std::vector<double>* TrkKeta = 0;
	std::vector<double>* TrkKphi = 0;	

	std::vector<double>* Trkpipt = 0;
	std::vector<double>* Trkpinhits = 0;
	std::vector<double>* Trkpichi2 = 0;
	std::vector<double>* Trkpidxy = 0;
	std::vector<double>* Trkpidz = 0;
	std::vector<double>* Trkpieta = 0;
	std::vector<double>* Trkpiphi = 0;

	std::vector<double>* TrkSpt = 0;
	std::vector<double>* TrkSnhits = 0;
	std::vector<double>* TrkSchi2 = 0;
	std::vector<double>* TrkSdxy = 0;
	std::vector<double>* TrkSdz = 0;
	std::vector<double>* TrkSeta = 0;
	std::vector<double>* TrkSphi = 0;
	std::vector<double>* DSDeltaR = 0;
	//----------------------------------------
	//For D* Wrong combination
	//-----------------------------------------
	std::vector<double>* D0massWrong = 0;
	std::vector<double>* D0etaWrong = 0;
	std::vector<double>* D0phiWrong = 0;
	std::vector<double>* D0ptWrong = 0;

	std::vector<double>* DsmassWrong = 0;
	std::vector<double>* DsetaWrong = 0;
	std::vector<double>* DsphiWrong = 0;
	std::vector<double>* DsptWrong = 0;

	//std::vector<double>* TrkScharge = 0;
	std::vector<double>* D0fromDSsXYWrong = 0;
	std::vector<double>* D0fromDSs3DWrong = 0;
	std::vector<double>* AnglephiWrong = 0;
	std::vector<double>* D0_VtxProbWrong = 0;
	//----------------------------------------
	//For D0 Closeth
	//-----------------------------------------
	std::vector<double>* D0Kpimass = 0;
	std::vector<double>* D0Kpi_VtxProb = 0;
	std::vector<double>* D0Kpis3D = 0;
	std::vector<double>* D0KpiDispAngle = 0;

	//=============================================
	//Variables to fill the new tree
	//=============================================
	string Variable_NameOfFiredTriggers; Variable_NameOfFiredTriggers.clear();
	
	//----------------------------------------
	//For D* 
	//-----------------------------------------
	double Variable_D0mass = 0.;
	double Variable_D0eta = 0.;
	double Variable_D0phi = 0.;
	double Variable_D0pt = 0.;
	double Variable_DsMinusD0 = 0.;
	double Variable_Dsmass = 0.;
	double Variable_Dseta = 0.;
	double Variable_Dsphi = 0.;
	double Variable_Dspt = 0.;

	double Variable_DSDeltaR = 0.;
	//double Variable_TrkScharge = 0.;
	double Variable_D0fromDSsXY = 0.;
	double Variable_D0fromDSs3D = 0.;
	double Variable_Anglephi = 0.;
	double Variable_D0_VtxProb = 0.;

	double Variable_TrkKpt = 0.;
	double Variable_TrkKnhits = 0.;
	double Variable_TrkKchi2 = 0.;
	double Variable_TrkKdxy = 0.;
	double Variable_TrkKdz = 0.;
	double Variable_TrkKeta = 0.;
	double Variable_TrkKphi = 0.;

	double Variable_Trkpipt = 0.;
	double Variable_Trkpinhits = 0.;
	double Variable_Trkpichi2 = 0.;
	double Variable_Trkpidxy = 0.;
	double Variable_Trkpidz = 0.;
	double Variable_Trkpieta = 0.;
	double Variable_Trkpiphi = 0.;

	double Variable_TrkSpt = 0.;
	double Variable_TrkSnhits = 0.;
	double Variable_TrkSchi2 = 0.;
	double Variable_TrkSdxy = 0.;
	double Variable_TrkSdz = 0.;
	double Variable_TrkSeta = 0.;
	double Variable_TrkSphi = 0.;
	//----------------------------------------
	//For D* Wrong combination
	//-----------------------------------------
	double Variable_D0massWrong = 0.;
	double Variable_D0etaWrong = 0.;
	double Variable_D0phiWrong = 0.;
	double Variable_D0ptWrong = 0.;
	double Variable_DsMinusD0Wrong = 0.;
	double Variable_DsmassWrong = 0.;
	double Variable_DsetaWrong = 0.;
	double Variable_DsphiWrong = 0.;
	double Variable_DsptWrong = 0.;

	//double Variable_TrkSchargeWrong = 0.;
	double Variable_D0fromDSsXYWrong = 0.;
	double Variable_D0fromDSs3DWrong = 0.;
	double Variable_AnglephiWrong = 0.;
	double Variable_D0_VtxProbWrong = 0.;

	//----------------------------------------
	//For D0 Closeth
	//-----------------------------------------
	double Variable_D0Kpimass = 0.;
	double Variable_D0Kpi_VtxProb = 0.;
	double Variable_D0Kpis3D = 0.;
	double Variable_D0KpiDispAngle = 0.;

	//VARIABLES D*

	t_analysis.Branch("NameOfFiredTriggers",&Variable_NameOfFiredTriggers, "Variable_NameOfFiredTriggers/S");
	t_analysis.Branch("D0mass",&Variable_D0mass, "Variable_D0mass/D");
	t_analysis.Branch("D0pt",&Variable_D0pt, "Variable_D0pt/D");
	t_analysis.Branch("D0eta",&Variable_D0eta, "Variable_D0eta/D");
	t_analysis.Branch("D0phi",&Variable_D0phi, "Variable_D0phi/D");
	t_analysis.Branch("DsMinusD0",&Variable_DsMinusD0, "Variable_DsMinusD0/D");
	t_analysis.Branch("Dsmass",&Variable_Dsmass, "Variable_Dsmass/D");
	t_analysis.Branch("Dspt",&Variable_Dspt, "Variable_Dspt/D");
	t_analysis.Branch("Dseta",&Variable_Dseta, "Variable_Dseta/D");
	t_analysis.Branch("Dsphi",&Variable_Dsphi, "Variable_Dsphi/D");

	t_analysis.Branch("DSDeltaR",&Variable_DSDeltaR, "Variable_DSDeltaR/D");
	t_analysis.Branch("D0_VtxProb",&Variable_D0_VtxProb, "Variable_D0_VtxProb/D");
	t_analysis.Branch("Anglephi",&Variable_Anglephi, "Variable_Anglephi/D");
	t_analysis.Branch("D0fromDSsXY",&Variable_D0fromDSsXY, "Variable_D0fromDSsXY/D");
	t_analysis.Branch("D0fromDSs3D",&Variable_D0fromDSs3D, "Variable_D0fromDSs3D/D");

	t_analysis.Branch("TrkKpt",&Variable_TrkKpt, "Variable_TrkKpt/D");
	t_analysis.Branch("Trkpipt",&Variable_Trkpipt, "Variable_Trkpipt/D");
	t_analysis.Branch("TrkSpt",&Variable_TrkSpt, "Variable_TrkSpt/D");
	t_analysis.Branch("TrkKnhits",&Variable_TrkKnhits, "Variable_TrkKnhits/D");
	t_analysis.Branch("Trkpinhits",&Variable_Trkpinhits, "Variable_Trkpinhits/D");
	t_analysis.Branch("TrkSnhits",&Variable_TrkSnhits, "Variable_TrkSnhits/D");
	t_analysis.Branch("TrkKchi2",&Variable_TrkKchi2, "Variable_TrkKchi2/D");
	t_analysis.Branch("Trkpichi2",&Variable_Trkpichi2, "Variable_Trkpichi2/D");
	t_analysis.Branch("TrkSchi2",&Variable_TrkSchi2, "Variable_TrkSchi2/D");
	t_analysis.Branch("TrkKdxy",&Variable_TrkKdxy, "Variable_TrkKdxy/D");
	t_analysis.Branch("Trkpidxy",&Variable_Trkpidxy, "Variable_Trkpidxy/D");
	t_analysis.Branch("TrkSdxy",&Variable_TrkSdxy, "Variable_TrkSdxy/D");
	t_analysis.Branch("TrkKdz",&Variable_TrkKdz, "Variable_TrkKdz/D");
	t_analysis.Branch("Trkpidz",&Variable_Trkpidz, "Variable_Trkpidz/D");
	t_analysis.Branch("TrkSdz",&Variable_TrkSdz, "Variable_TrkSdz/D"); 
	t_analysis.Branch("TrkKeta",&Variable_TrkKeta, "Variable_TrkKeta/D");      
	t_analysis.Branch("Trkpieta",&Variable_Trkpieta, "Variable_Trkpieta/D");
	t_analysis.Branch("TrkSeta",&Variable_TrkSeta, "Variable_TrkSeta/D");
	t_analysis.Branch("TrkKphi",&Variable_TrkKphi, "Variable_TrkKphi/D");
	t_analysis.Branch("Trkpiphi",&Variable_Trkpiphi, "Variable_Trkpiphi/D");
	t_analysis.Branch("TrkSphi",&Variable_TrkSphi, "Variable_TrkSphi/D");
	//t_analysis.Branch("TrkScharge",&Variable_TrkScharge, "Variable_TrkScharge/D");

	//----------------------------------------
	//For D* Wrong combination
	//-----------------------------------------
	t_DstarWrongCombination.Branch("D0massWrong",&Variable_D0massWrong, "Variable_D0massWrong/D");
	t_DstarWrongCombination.Branch("D0ptWrong",&Variable_D0ptWrong, "Variable_D0ptWrong/D");
	t_DstarWrongCombination.Branch("D0etaWrong",&Variable_D0etaWrong, "Variable_D0etaWrong/D");
	t_DstarWrongCombination.Branch("D0phiWrong",&Variable_D0phiWrong, "Variable_D0phiWrong/D");

	t_DstarWrongCombination.Branch("DsMinusD0Wrong",&Variable_DsMinusD0Wrong, "Variable_DsMinusD0Wrong/D");
	
	t_DstarWrongCombination.Branch("DsmassWrong",&Variable_DsmassWrong, "Variable_DsmassWrong/D");
	t_DstarWrongCombination.Branch("DsptWrong",&Variable_DsptWrong, "Variable_DsptWrong/D");
	t_DstarWrongCombination.Branch("DsetaWrong",&Variable_DsetaWrong, "Variable_DsetaWrong/D");
	t_DstarWrongCombination.Branch("DsphiWrong",&Variable_DsphiWrong, "Variable_DsphiWrong/D");

	t_DstarWrongCombination.Branch("D0_VtxProbWrong",&Variable_D0_VtxProbWrong, "Variable_D0_VtxProbWrong/D");
	t_DstarWrongCombination.Branch("AnglephiWrong",&Variable_AnglephiWrong, "Variable_AnglephiWrong/D");
	t_DstarWrongCombination.Branch("D0fromDSsXYWrong",&Variable_D0fromDSsXYWrong, "Variable_D0fromDSsXYWrong/D");
	t_DstarWrongCombination.Branch("D0fromDSs3DWrong",&Variable_D0fromDSs3DWrong, "Variable_D0fromDSs3DWrong/D");
	
	//----------------------------------------
	//For D0 Closeth
	//-----------------------------------------
	t_D0analysis.Branch("D0Kpimass",&Variable_D0Kpimass, "Variable_D0Kpimass/D");
	t_D0analysis.Branch("D0Kpi_VtxProb",&Variable_D0Kpi_VtxProb, "Variable_D0Kpi_VtxProb/D");
	t_D0analysis.Branch("D0Kpis3D",&Variable_D0Kpis3D, "Variable_D0Kpis3D/D");
	t_D0analysis.Branch("D0KpiDispAngle",&Variable_D0KpiDispAngle, "Variable_D0KpiDispAngle/D");

	if (debug)cout << "debug 2 --------------------" << endl;
	//--------------------------------------------------
	//Creating Histgrams
	//---------------------------------------------------
	//TH1* makeTH1(const char* name, Double_t BIN, Double_t NMIN, Double_t NMAX, const char* TITLES, Color_t COLOR,const char* TYPE)
	TH1 *D0massHisto = makeTH1("D0massHisto", 100, 1.76, 1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassHisto = makeTH1("DsmassHisto", 100,1.93,2.12, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *D0fromDSs3DHisto = makeTH1("D0fromDSs3DHisto", 100,0,5, "Significance of the D0 ; Significance ; Events ", kRed);
	TH1 *TotalD0fromDSs3DHisto = makeTH1("TotalD0fromDSs3DHisto", 100,0,5, "Significance of the D0 ; Significance ; Events ", kRed);
	TH1 *DsMinusD0Histo = makeTH1("DsMinusD0Histo", 100,0.14,0.16, "#Deltam = m(K#pi#pi_{slow}) - m(K#pi) ; #Delta m ; Events", kRed);
	TH1 *D0VtxProbHisto = makeTH1("D0VtxProbHisto", 100,0,1, "D0VtxProbHisto ; Probability ; Events ", kRed);

	TH1 *D0KpimassHisto = makeTH1("D0KpimassHisto", 100,1.76,1.96, "Invariant Mass of the D0(Prompt) ; Mass [GeV] ; Events ", kRed);
	TH1 *D0Kpi_VtxProbHisto = makeTH1("D0Kpi_VtxProbHisto", 100,0,1, "D0Kpi_VtxProbHisto D0(Prompt) ; Probability ; Events ", kRed);
	TH1 *D0Kpis3DHisto = makeTH1("D0Kpis3DHisto", 100,0,5, "Significance of the D0(Prompt) ; Significance ; Events ", kRed);
	//TH1 *D0KpiDispAngleHisto = makeTH1("D0KpiDispAngleHisto", 100,0,1, "Invariant Mass of the D0(Prompt) ; Mass [GeV] ; Events ", kRed);

	if (debug) cout << "debug 3 --------------------"<< endl;
	
	//--------------------------------------------------
	//ADDRESSING THE MEMORY TO VECTOR AND VARIABLES
	//--------------------------------------------------
	TBranch *b_NameOfFiredTriggers; t1->SetBranchAddress("NameOfFiredTriggers",&NameOfFiredTriggers,&b_NameOfFiredTriggers);
	//For D* 
	//-----------------------------------------
	//TBranch *b_D0mass = t1->GetBranch("D0mass");b_D0mass->SetAddress(&D0mass);
	TBranch *b_D0mass; t1->SetBranchAddress("D0mass",&D0mass,&b_D0mass);
	TBranch *b_Dsmass; t1->SetBranchAddress("Dsmass",&Dsmass,&b_Dsmass);
	TBranch *b_D0eta; t1->SetBranchAddress("D0eta",&D0eta,&b_D0eta);
	TBranch *b_D0phi; t1->SetBranchAddress("D0phi",&D0phi,&b_D0phi);
	TBranch *b_D0pt; t1->SetBranchAddress("D0pt",&D0pt,&b_D0pt);

	TBranch *b_Dseta; t1->SetBranchAddress("Dseta",&Dseta,&b_Dseta);
	TBranch *b_Dsphi ; t1->SetBranchAddress("Dsphi",&Dsphi,&b_Dsphi);
	TBranch *b_Dspt ; t1->SetBranchAddress("Dspt",&Dspt,&b_Dspt);
	TBranch *b_DSDeltaR ; t1->SetBranchAddress("DSDeltaR",&DSDeltaR,&b_DSDeltaR);
	TBranch *b_Anglephi ; t1->SetBranchAddress("Anglephi",&Anglephi,&b_Anglephi);
	TBranch *b_D0fromDSsXY ; t1->SetBranchAddress("D0fromDSsXY",&D0fromDSsXY,&b_D0fromDSsXY);
	TBranch *b_D0fromDSs3D ; t1->SetBranchAddress("D0fromDSs3D",&D0fromDSs3D,&b_D0fromDSs3D);
	TBranch *b_D0_VtxProb ; t1->SetBranchAddress("D0_VtxProb",&D0_VtxProb,&b_D0_VtxProb);

	TBranch *b_TrkSpt ; t1->SetBranchAddress("TrkSpt",&TrkSpt,&b_TrkSpt);
	TBranch *b_TrkSnhits ; t1->SetBranchAddress("TrkSnhits",&TrkSnhits,&b_TrkSnhits);
	TBranch *b_TrkSchi2 ; t1->SetBranchAddress("TrkSchi2",&TrkSchi2,&b_TrkSchi2);
	TBranch *b_TrkSdxy ; t1->SetBranchAddress("TrkSdxy",&TrkSdxy,&b_TrkSdxy);
	TBranch *b_TrkSdz ; t1->SetBranchAddress("TrkSdz",&TrkSdz,&b_TrkSdz);
	TBranch *b_TrkSeta ; t1->SetBranchAddress("TrkSeta",&TrkSeta,&b_TrkSeta);
	TBranch *b_TrkSphi ; t1->SetBranchAddress("TrkSphi",&TrkSphi,&b_TrkSphi);

	TBranch *b_Trkpipt ; t1->SetBranchAddress("Trkpipt",&Trkpipt,&b_Trkpipt);
	TBranch *b_Trkpinhits ; t1->SetBranchAddress("Trkpinhits",&Trkpinhits,&b_Trkpinhits);
	TBranch *b_Trkpichi2 ; t1->SetBranchAddress("Trkpichi2",&Trkpichi2,&b_Trkpichi2);
	TBranch *b_Trkpidxy ; t1->SetBranchAddress("Trkpidxy",&Trkpidxy,&b_Trkpidxy);
	TBranch *b_Trkpidz ; t1->SetBranchAddress("Trkpidz",&Trkpidz,&b_Trkpidz);
	TBranch *b_Trkpieta ; t1->SetBranchAddress("Trkpieta",&Trkpieta,&b_Trkpieta);
	TBranch *b_Trkpiphi ; t1->SetBranchAddress("Trkpiphi",&Trkpiphi,&b_Trkpiphi);

	TBranch *b_TrkKpt ; t1->SetBranchAddress("TrkKpt",&TrkKpt,&b_TrkKpt);
	TBranch *b_TrkKnhits ; t1->SetBranchAddress("TrkKnhits",&TrkKnhits,&b_TrkKnhits);
	TBranch *b_TrkKchi2 ; t1->SetBranchAddress("TrkKchi2",&TrkKchi2,&b_TrkKchi2);
	TBranch *b_TrkKdxy ; t1->SetBranchAddress("TrkKdxy",&TrkKdxy,&b_TrkKdxy);
	TBranch *b_TrkKdz ; t1->SetBranchAddress("TrkKdz",&TrkKdz,&b_TrkKdz);
	TBranch *b_TrkKeta ; t1->SetBranchAddress("TrkKeta",&TrkKeta,&b_TrkKeta);
	TBranch *b_TrkKphi ; t1->SetBranchAddress("TrkKphi",&TrkKphi,&b_TrkKphi);

	//----------------------------------------
	//For D* WRONG COMBINATION
	//-----------------------------------------
	TBranch *b_D0massWrong; t1->SetBranchAddress("D0massWrong",&D0massWrong,&b_D0massWrong);
	TBranch *b_DsmassWrong; t1->SetBranchAddress("DsmassWrong",&DsmassWrong,&b_DsmassWrong);
	TBranch *b_D0etaWrong; t1->SetBranchAddress("D0etaWrong",&D0etaWrong,&b_D0etaWrong);
	TBranch *b_D0phiWrong; t1->SetBranchAddress("D0phiWrong",&D0phiWrong,&b_D0phiWrong);
	TBranch *b_D0ptWrong; t1->SetBranchAddress("D0ptWrong",&D0ptWrong,&b_D0ptWrong);

	TBranch *b_DsetaWrong ; t1->SetBranchAddress("DsetaWrong",&DsetaWrong,&b_DsetaWrong);
	TBranch *b_DsphiWrong ; t1->SetBranchAddress("DsphiWrong",&DsphiWrong,&b_DsphiWrong);
	TBranch *b_DsptWrong ; t1->SetBranchAddress("DsptWrong",&DsptWrong,&b_DsptWrong);
	TBranch *b_AnglephiWrong ; t1->SetBranchAddress("AnglephiWrong",&AnglephiWrong,&b_AnglephiWrong);
	TBranch *b_D0fromDSsXYWrong ; t1->SetBranchAddress("D0fromDSsXYWrong",&D0fromDSsXYWrong,&b_D0fromDSsXYWrong);
	TBranch *b_D0fromDSs3DWrong ; t1->SetBranchAddress("D0fromDSs3DWrong",&D0fromDSs3DWrong,&b_D0fromDSs3DWrong);
	TBranch *b_D0_VtxProbWrong ; t1->SetBranchAddress("D0_VtxProbWrong",&D0_VtxProbWrong,&b_D0_VtxProbWrong);

	//----------------------------------------
	//FOR D0 CLOSETH
	//-----------------------------------------
	TBranch *b_D0Kpimass ; t1->SetBranchAddress("D0Kpimass",&D0Kpimass,&b_D0Kpimass);
	TBranch *b_D0Kpi_VtxProb ; t1->SetBranchAddress("D0Kpi_VtxProb",&D0Kpi_VtxProb,&b_D0Kpi_VtxProb);
	TBranch *b_D0Kpis3D ; t1->SetBranchAddress("D0Kpis3D",&D0Kpis3D,&b_D0Kpis3D);
	TBranch *b_D0KpiDispAngle ; t1->SetBranchAddress("D0KpiDispAngle",&D0KpiDispAngle,&b_D0KpiDispAngle);
	//TBranch *b_teste = t1->GetBranch("teste");b_teste->SetAddress(&teste);	
	if (debug)cout << "debug 5 --------------------"<< endl;

	//**********************************************************		
	//Reading Number of tree entries of the file
	Long64_t nentries = t1->GetEntries();
	cout<< "Number of tree entries: "<< nentries <<std::endl;
	//Long64_t GetEntriesFast = t1->GetEntriesFast();

	for (Long64_t jentry=0; jentry < nentries; jentry++) //loop tree entries for file f1
	{
		if (debug)cout << "debug 7 --------------------"<< endl;
		Long64_t ientry = t1->LoadTree(jentry);
      //std::cout << "nentries " << nentries << " ientry " << ientry << " jentry " << jentry <<std::endl;
      if (ientry < 0) break;
		
		//Output about percent program executed
		double percent = (jentry*100)/nentries;
		if (jentry % 10000 == 0) cout <<"*****"<< percent << "per cent done***" << endl;
		//cout <<"*****entrada"<< jentry << "***" << endl;
		

		b_NameOfFiredTriggers->GetEntry(ientry);
		//----------------------------------------
		//For D* 
		//-----------------------------------------
		b_D0mass->GetEntry(ientry);
		b_D0eta->GetEntry(ientry);
		b_D0phi->GetEntry(ientry);
		b_D0pt->GetEntry(ientry);
		b_Dsmass->GetEntry(ientry);
		b_Dseta->GetEntry(ientry);
		b_Dsphi->GetEntry(ientry);
		b_Dspt->GetEntry(ientry);
		b_Anglephi->GetEntry(ientry);
		b_D0fromDSsXY->GetEntry(ientry);
		b_D0fromDSs3D->GetEntry(ientry);
		b_D0_VtxProb->GetEntry(ientry);
		b_DSDeltaR->GetEntry(ientry);

		b_TrkSpt->GetEntry(ientry);
		b_TrkSnhits->GetEntry(ientry);
		b_TrkSchi2->GetEntry(ientry);
		b_TrkSdxy->GetEntry(ientry);
		b_TrkSdz->GetEntry(ientry);
		b_TrkSeta->GetEntry(ientry);
		b_TrkSphi->GetEntry(ientry);

		b_TrkKpt->GetEntry(ientry);
		b_TrkKnhits->GetEntry(ientry);
		b_TrkKchi2->GetEntry(ientry);
		b_TrkKdxy->GetEntry(ientry);
		b_TrkKdz->GetEntry(ientry);
		b_TrkKeta->GetEntry(ientry);
		b_TrkKphi->GetEntry(ientry);

		b_Trkpipt->GetEntry(ientry);
		b_Trkpinhits->GetEntry(ientry);
		b_Trkpichi2->GetEntry(ientry);
		b_Trkpidxy->GetEntry(ientry);
		b_Trkpidz->GetEntry(ientry);
		b_Trkpieta->GetEntry(ientry);
		b_Trkpiphi->GetEntry(ientry);
		//----------------------------------------
		//For D* Wrong combination
		//-----------------------------------------
		b_D0massWrong->GetEntry(ientry);
		b_D0etaWrong->GetEntry(ientry);
		b_D0phiWrong->GetEntry(ientry);
		b_D0ptWrong->GetEntry(ientry);	
		b_DsmassWrong->GetEntry(ientry);
		b_DsetaWrong->GetEntry(ientry);
		b_DsphiWrong->GetEntry(ientry);
		b_DsptWrong->GetEntry(ientry);
		b_AnglephiWrong->GetEntry(ientry);
		b_D0fromDSsXYWrong->GetEntry(ientry);
		b_D0fromDSs3DWrong->GetEntry(ientry);
		b_D0_VtxProbWrong->GetEntry(ientry);
		//----------------------------------------
		//For D0 Closeth
		//-----------------------------------------
		b_D0Kpimass->GetEntry(ientry);
		b_D0Kpi_VtxProb->GetEntry(ientry);
		b_D0Kpis3D->GetEntry(ientry);
		b_D0KpiDispAngle->GetEntry(ientry);
	
		if (debug)cout << "debug 10 --------------------"<< endl;

		//if (std::find(NameOfFiredTriggers->begin(), NameOfFiredTriggers->end(), "HLT_Mu9_IP6_part0") == NameOfFiredTriggers->end()) continue;
		//if (std::find_if(NameOfFiredTriggers->begin(), NameOfFiredTriggers->end(), NameOfFiredTriggers.find("HLT_Mu9_IP6") ) continue;
		std::cout << "--------------------"<< std::endl;

		for(unsigned int tescont2=0; tescont2 < NameOfFiredTriggers->size(); tescont2++)
		{  		
			string teststring = NameOfFiredTriggers->at(tescont);
			std::cout << "NameOfFiredTriggers ["<< tescont <<"]: " << teststring  << std::endl;
		}

		//----------------------------------------
		//For D* and D0(from D*)
		//-----------------------------------------	
		for(unsigned int i=0; i < D0mass->size(); i++)
		{  if (debug)cout << "debug D* 11 --------------------"<< endl;
							
			if( D0_VtxProb->at(i) < 0.01) continue;
			TotalD0fromDSs3DHisto->Fill(D0fromDSs3D->at(i));
			if( D0fromDSs3D->at(i) < SigCutD0fromDs) continue;
			if( Anglephi->at(i) < 0.99 ) continue;
			if ( D0pt->at(i) < 3. ) continue;
			if( (Dsmass->at(i) - D0mass->at(i)) > 0.16) continue;
						
			D0massHisto->Fill(D0mass->at(i));
			DsmassHisto->Fill(Dsmass->at(i));
			DsMinusD0Histo->Fill((Dsmass->at(i) - D0mass->at(i)));
			D0fromDSs3DHisto->Fill(D0fromDSs3D->at(i));
			D0VtxProbHisto->Fill(D0_VtxProb->at(i));

			Variable_DsMinusD0 = (Dsmass->at(i) - D0mass->at(i));
	
			Variable_D0mass = D0mass->at(i);
			Variable_Dsmass = Dsmass->at(i);
			Variable_D0eta = D0eta->at(i);
			Variable_D0phi = D0phi->at(i);
			Variable_Dseta = Dseta->at(i);
			Variable_Dsphi = Dsphi->at(i);
			Variable_TrkKpt = TrkKpt->at(i);
			Variable_D0pt = D0pt->at(i);
			Variable_Dspt = Dspt->at(i);
			Variable_Trkpipt = Trkpipt->at(i);
			Variable_TrkSpt = TrkSpt->at(i);
			Variable_DSDeltaR = DSDeltaR->at(i);
			Variable_TrkKnhits = TrkKnhits->at(i);
			Variable_Trkpinhits = Trkpinhits->at(i);
			Variable_TrkSnhits = TrkSnhits->at(i);
			Variable_TrkKchi2 = TrkKchi2->at(i);
			Variable_Trkpichi2 = Trkpichi2->at(i);
			Variable_TrkSchi2 = TrkSchi2->at(i);
			Variable_TrkKdxy = TrkKdxy->at(i);
			Variable_Trkpidxy = Trkpidxy->at(i);
			Variable_TrkSdxy = TrkSdxy->at(i);
			Variable_TrkKdz = TrkKdz->at(i);
			Variable_Trkpidz = Trkpidz->at(i);
			Variable_TrkSdz = TrkSdz->at(i);
			Variable_TrkKeta = TrkKeta->at(i);
			Variable_Trkpieta = Trkpieta->at(i);
			Variable_TrkSeta = TrkSeta->at(i);
			Variable_TrkKphi = TrkKphi->at(i);
			Variable_Trkpiphi = Trkpiphi->at(i);
			Variable_TrkSphi = TrkSphi->at(i);
			//Variable_TrkScharge = TrkScharge->at(i);
			Variable_D0fromDSsXY = D0fromDSsXY->at(i);
			Variable_D0fromDSs3D = D0fromDSs3D->at(i);
			Variable_Anglephi = Anglephi->at(i);
			Variable_D0_VtxProb = D0_VtxProb->at(i);

			t_analysis.Fill();

			if (debug)cout << "debug D* 12 --------------------"<< endl;
  		}
	
		//----------------------------------------
		//For D* Wrong Combination
		//-----------------------------------------	
		for(unsigned int k=0; k < D0massWrong->size(); k++)
		{ 				
			if( D0_VtxProbWrong->at(k) < 0.01) continue;
			if( D0fromDSs3DWrong->at(k) < SigCutD0fromDs) continue;
			if( AnglephiWrong->at(k) < 0.99 ) continue;
			if ( D0ptWrong->at(k) < 3. ) continue;
			if( (DsmassWrong->at(k) - D0massWrong->at(k)) > 0.16) continue;
						
			Variable_DsMinusD0 = (DsmassWrong->at(k) - D0massWrong->at(k));
			Variable_D0massWrong = D0massWrong->at(k);
			Variable_DsmassWrong = DsmassWrong->at(k);
			Variable_D0etaWrong = D0etaWrong->at(k);
			Variable_D0phiWrong = D0phiWrong->at(k);
			Variable_DsetaWrong = DsetaWrong->at(k);
			Variable_DsphiWrong = DsphiWrong->at(k);
			Variable_D0ptWrong = D0ptWrong->at(k);
			Variable_DsptWrong = DsptWrong->at(k);

			t_DstarWrongCombination.Fill();	
  		}//For D* Wrong Combination

		//----------------------------------------
		//For D0(Direct)
		//-----------------------------------------	
		for(unsigned int j=0; j < D0Kpimass->size(); j++)
		{  if (debug)cout << "debug D0(Direct) 11 --------------------"<< endl;
							
			if( D0Kpi_VtxProb->at(j) < 0.01) continue;
			if( D0Kpis3D->at(j) < SigCutD0) continue;
			if( D0KpiDispAngle->at(j) < 0.99 ) continue;
			//if ( D0Kpipt->at(i) < 3. ) continue;	
		
			D0KpimassHisto->Fill(D0Kpimass->at(j));
			D0Kpis3DHisto->Fill(D0Kpis3D->at(j));
			D0Kpi_VtxProbHisto->Fill(D0Kpi_VtxProb->at(j));

			Variable_D0Kpimass = D0Kpimass->at(j);
			Variable_D0Kpi_VtxProb = D0Kpi_VtxProb->at(j);
			Variable_D0Kpis3D = D0Kpis3D->at(j);
			Variable_D0KpiDispAngle = D0KpiDispAngle->at(j);
		
			//vectorInvariantMass_D0.push_back(D0mass->at(i));
			if (debug)cout << "debug D0(Direct) 12 --------------------"<< endl;
			t_D0analysis.Fill();	
  		}//For D0(Direct)
			
		if (debug)cout << "debug 13 --------------------"<< endl;	
	}//End loop tree entries for file f1

	if (debug)cout << "debug 14 --------------------"<< endl;

	TCanvas* canvas = new TCanvas("canvas","",900,600);
	canvas->Divide(3,3);
	canvas->cd(1);
	D0massHisto->Draw("e1p");
	canvas->cd(2);
	DsmassHisto->Draw("e1p");
	canvas->cd(3);
	//TotalD0fromDSs3DHisto->Draw("Histo");
	D0fromDSs3DHisto->Draw("Histo");
	canvas->cd(4);	
	DsMinusD0Histo->Draw("Histo");
	canvas->cd(5);
	D0VtxProbHisto->Draw("Histo");
	canvas->cd(6);
	D0KpimassHisto->Draw("e1p");
	canvas->cd(7);
	D0Kpis3DHisto->Draw("Histo");
	canvas->cd(8);
	D0Kpi_VtxProbHisto->Draw("Histo");

	canvas->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/PreliminarStudy.pdf");

	if (debug)cout << "debug 15 --------------------"<< endl;

	//t_analysis.Branch("D0mass",&D0mass);

	//TBranch *D0mass_branch; D0mass_branch = 
	f_analysis.cd();
	D0massHisto->Write();//Write() if a file is open, this function writes a root objectics on it.
	DsmassHisto->Write();
	DsMinusD0Histo->Write();
	t_analysis.Write();  //Write in the root file

	t_DstarWrongCombination.Write();

	D0KpimassHisto->Write();
	t_D0analysis.Write();

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout << "-------END PROGRAM-------------"<< endl;

}//end program
