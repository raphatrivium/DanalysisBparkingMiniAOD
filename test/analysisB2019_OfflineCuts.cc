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
//#include "Riosfwd.h"
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

	string Dataset = "MC_DStarToD0Pi_D0KPi"; 
	//Datasets: BparkingDataset1 - BparkingDataset2 - BparkingDataset3
	//Datasets: BparkingDataset4 - BparkingDataset5 - BparkingDataset6
	//Datasets: MC_DStarToD0Pi_D0KPi

	string TriggerPath = ""; //HLT_Mu9_IP6
	 
	//Save the file
	string path = "/eos/user/r/ragomesd/crab/haddteste/";
	double SigCutD0fromDs = 1.0; //3.
	double SigCutD0 = 1.0; //3.

	bool debug = false;
	//eos/user/r/ragomesd/crab
		
	//TChain* t1 = callTchain("/eos/user/r/ragomesd/crab/DStarToD0Pi_D0KPi_NoMuFilter_TuneCP5_13TeV-pythia8-evtgen/Ds_MC_Run2018A_Official/191215_044633/0000/");
	//TChain* t1 = callOneFile("/eos/user/r/ragomesd/testes/D0DstarDataBparking.root");
	//TChain* t1 = new TChain("analysis/data"); t1->Add("/eos/user/r/ragomesd/crab/ParkingBPH1/ParkingBPH_Run2018A_MINIAOD/191212_023520/0000/D0DstarDataBparking_35.root");
	TChain* t1 = new TChain("analysis/data"); t1->Add("/afs/cern.ch/work/r/ragomesd/B_parking/CMSSW_10_2_7/src/BanalysisMiniAOD-1/DanalysisBparkingMiniAOD/python/D0DstarMC.root");
	
	
	//-------Reading the root file-------------------------------------	
	//TFile *f1 = new TFile("D0DstarDataBparking_7.root");
	//TTree *t1 = (TTree*)f1->Get("analysis/data");
	
	//-------Create root file-------------------------------------
	//TFile f_analysis("DMesonsHistograms.root","recreate"); //Creates root file
	TFile f_analysis(Form("%s%s.root", path.c_str(), Dataset.c_str() ),"recreate"); //Creates root file
	TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	TTree t_DstarWrongCombination("t_DstarWrongCombination","t_DstarWrongCombination");
	TTree t_D0analysis("t_D0analysis","t_D0analysis");
	TTree t_DsMC("t_DsMC","t_DsMC"); 
	TTree t_D0MC("t_D0MC","t_D0MC"); 
	TTree t_DsMatching("t_DsMatching","t_DsMatching");
	TTree t_D0Matching("t_D0Matching","t_D0Matching"); 

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
	//For D0
	//-----------------------------------------
	std::vector<double>* D0Kpimass = 0;
	std::vector<double>* D0Kpi_VtxProb = 0;
	std::vector<double>* D0Kpis3D = 0;
	std::vector<double>* D0KpiDispAngle = 0;

	std::vector<double>* D0Kpipt = 0;
	std::vector<double>* D0Kpieta = 0;
	std::vector<double>* D0Kpiphi = 0;

	std::vector<double>* TrkD0Kdxy = 0;
	std::vector<double>* TrkD0Kdz = 0;
	std::vector<double>* TrkD0Kchi2 = 0;
	std::vector<double>* TrkD0Kpt = 0;
	std::vector<double>* TrkD0Keta = 0;
	std::vector<double>* TrkD0kphi = 0;
	std::vector<double>* TrkD0Knhits = 0;

	std::vector<double>* TrkD0pidxy = 0;
	std::vector<double>* TrkD0pidz = 0;
	std::vector<double>* TrkD0pichi2 = 0;
	std::vector<double>* TrkD0pipt = 0;
	std::vector<double>* TrkD0pieta = 0;
	std::vector<double>* TrkD0piphi = 0;
	std::vector<double>* TrkD0pinhits = 0;

	std::vector<double>* D0KpisXY = 0;
	std::vector<double>* D0Kpid3D = 0;
	std::vector<double>* D0Kpie3D = 0;
	std::vector<double>* D0_kT = 0;

	double Variable_D0Kpimass = 0.;
	double Variable_D0Kpi_VtxProb = 0.;
	double Variable_D0Kpis3D = 0.;
	double Variable_D0KpiDispAngle = 0.;

	double Variable_D0Kpipt = 0.;
	double Variable_D0Kpieta = 0.;
	double Variable_D0Kpiphi = 0.;

	double Variable_TrkD0Kdxy = 0.;
	double Variable_TrkD0Kdz = 0.;
	double Variable_TrkD0Kchi2 = 0.;
	double Variable_TrkD0Kpt = 0.;
	double Variable_TrkD0Keta = 0.;
	double Variable_TrkD0kphi = 0.;
	double Variable_TrkD0Knhits = 0.;

	double Variable_TrkD0pidxy = 0.;
	double Variable_TrkD0pidz = 0.;
	double Variable_TrkD0pichi2 = 0.;
	double Variable_TrkD0pipt = 0.;
	double Variable_TrkD0pieta = 0.;
	double Variable_TrkD0piphi = 0.;
	double Variable_TrkD0pinhits = 0.;

	double Variable_D0KpisXY = 0.;
	double Variable_D0Kpid3D = 0.;
	double Variable_D0Kpie3D = 0.;
	double Variable_D0_kT = 0.;

//======================================================
// MC Variables - D_star
//======================================================
	std::vector<double>* dScandsKpi = 0;
	std::vector<double>* MCDseta = 0;
	std::vector<double>* MCDsphi = 0;
	std::vector<double>* MCDspt = 0;
	std::vector<double>* MCDsenergy = 0;
	std::vector<double>* MCDsp = 0;
	std::vector<double>* MCDset = 0;
	std::vector<double>* MCDsmass = 0;
	std::vector<double>* MCD0eta = 0;
	std::vector<double>* MCD0phi = 0;
	std::vector<double>* MCD0pt = 0;
	std::vector<double>* MCD0energy = 0;
	std::vector<double>* MCD0p = 0;
	std::vector<double>* MCD0et = 0;
	std::vector<double>* MCD0rapidity = 0;
	std::vector<double>* MCD0mass = 0;
	std::vector<double>* MCDsKphi = 0;
	std::vector<double>* MCDsKpt = 0;
	std::vector<double>* MCDsKenergy = 0;
	std::vector<double>* MCDsKp = 0;
	std::vector<double>* MCDsKet = 0;
	std::vector<double>* MCDsKrapidity = 0;
	std::vector<double>* MCDsKmass = 0;
	std::vector<double>* MCDsPieta = 0;
	std::vector<double>* MCDsPiphi = 0;
	std::vector<double>* MCDsPipt = 0;
	std::vector<double>* MCDsPienergy = 0;
	std::vector<double>* MCDsPip = 0;
	std::vector<double>* MCDsPiet = 0;
	std::vector<double>* MCDsPirapidity = 0;
	std::vector<double>* MCDsPimass = 0;

	double Variable_dScandsKpi = 0.;
	double Variable_MCDseta = 0.;
	double Variable_MCDsphi = 0.;
	double Variable_MCDspt = 0.;
	double Variable_MCDsenergy = 0.;
	double Variable_MCDsp = 0.;
	double Variable_MCDset = 0.;
	double Variable_MCDsmass = 0.;
	double Variable_MCD0eta = 0.;
	double Variable_MCD0phi = 0.;
	double Variable_MCD0pt = 0.;
	double Variable_MCD0energy = 0.;
	double Variable_MCD0p = 0.;
	double Variable_MCD0et = 0.;
	double Variable_MCD0rapidity = 0.;
	double Variable_MCD0mass = 0.;
	double Variable_MCDsKphi = 0.;
	double Variable_MCDsKpt = 0.;
	double Variable_MCDsKenergy = 0.;
	double Variable_MCDsKp = 0.;
	double Variable_MCDsKet = 0.;
	double Variable_MCDsKrapidity = 0.;
	double Variable_MCDsKmass = 0.;
	double Variable_MCDsPieta = 0.;
	double Variable_MCDsPiphi = 0.;
	double Variable_MCDsPipt = 0.;
	double Variable_MCDsPienergy = 0.;
	double Variable_MCDsPip = 0.;
	double Variable_MCDsPiet = 0.;
	double Variable_MCDsPirapidity = 0.;
	double Variable_MCDsPimass = 0.;

//======================================================
// MC Variables - D0 
//======================================================
	std::vector<double>* MCpromptD0eta = 0;
	std::vector<double>* MCpromptD0phi = 0;
	std::vector<double>* MCpromptD0pt = 0;
	std::vector<double>* MCpromptD0energy = 0;
	std::vector<double>* MCpromptD0p = 0;
	std::vector<double>* MCpromptD0et = 0;
	std::vector<double>* MCpromptD0rapidity = 0;
	std::vector<double>* MCpromptD0mass = 0;
	std::vector<double>* MCpromptD0_Keta = 0;
	std::vector<double>* MCpromptD0_Kphi = 0;
	std::vector<double>* MCpromptD0_Kpt = 0;
	std::vector<double>* MCpromptD0_Kenergy = 0;
	std::vector<double>* MCpromptD0_Kp = 0;
	std::vector<double>* MCpromptD0_Ket = 0;
	std::vector<double>* MCpromptD0_Krapidity = 0;
	std::vector<double>* MCpromptD0_Kmass = 0;
	std::vector<double>* MCpromptD0_Pieta = 0;
	std::vector<double>* MCpromptD0_Piphi = 0;
	std::vector<double>* MCpromptD0_Pipt = 0;
	std::vector<double>* MCpromptD0_Pienergy = 0;
	std::vector<double>* MCpromptD0_Pip = 0;
	std::vector<double>* MCpromptD0_Piet = 0;
	std::vector<double>* MCpromptD0_Pirapidity = 0;
	std::vector<double>* MCpromptD0_Pimass = 0;
	std::vector<double>* MCpromptD0_DispAngle = 0;

	double Variable_MCpromptD0eta = 0.;
	double Variable_MCpromptD0phi = 0.;
	double Variable_MCpromptD0pt = 0.;
	double Variable_MCpromptD0energy = 0.;
	double Variable_MCpromptD0p = 0.;
	double Variable_MCpromptD0et = 0.;
	double Variable_MCpromptD0rapidity = 0.;
	double Variable_MCpromptD0mass = 0.;
	double Variable_MCpromptD0_Keta = 0.;
	double Variable_MCpromptD0_Kphi = 0.;
	double Variable_MCpromptD0_Kpt = 0.;
	double Variable_MCpromptD0_Kenergy = 0.;
	double Variable_MCpromptD0_Kp = 0.;
	double Variable_MCpromptD0_Ket = 0.;
	double Variable_MCpromptD0_Krapidity = 0.;
	double Variable_MCpromptD0_Kmass = 0.;
	double Variable_MCpromptD0_Pieta = 0.;
	double Variable_MCpromptD0_Piphi = 0.;
	double Variable_MCpromptD0_Pipt = 0.;
	double Variable_MCpromptD0_Pienergy = 0.;
	double Variable_MCpromptD0_Pip = 0.;
	double Variable_MCpromptD0_Piet = 0.;
	double Variable_MCpromptD0_Pirapidity = 0.;
	double Variable_MCpromptD0_Pimass = 0.;
	double Variable_MCpromptD0_DispAngle = 0.;

//======================================================
// MC Matching - D* 
//======================================================

	std::vector<double>* DsetaMatching = 0;
	std::vector<double>* MCDsetaMatching = 0;
	std::vector<double>* DsphiMatching = 0;
	std::vector<double>* MCDsphiMatching = 0;
	std::vector<double>* DsptMatching = 0;
	std::vector<double>* MCDsptMatching = 0; 
	std::vector<double>* D0fromDsmass = 0;
	std::vector<double>* deltaRDsMatching = 0;

	double Variable_DsetaMatching = 0.;
	double Variable_MCDsetaMatching = 0.;
	double Variable_DsphiMatching = 0.;
	double Variable_MCDsphiMatching = 0.;
	double Variable_DsptMatching = 0.;
	double Variable_MCDsptMatching = 0.;
	double Variable_D0fromDsmass = 0.;
	double Variable_deltaRDsMatching = 0.;

//======================================================
// MC Matching - D0 
//======================================================

	std::vector<double>* D0etaMatching = 0;
	std::vector<double>* MCD0etaMatching = 0;
	std::vector<double>* D0phiMatching = 0;
	std::vector<double>* MCD0phiMatching = 0;
	std::vector<double>* D0ptMatching = 0;
	std::vector<double>* MCD0ptMatching = 0;
	std::vector<double>* D0KtMatching = 0;
	std::vector<double>* D0SxyMatching = 0;
	std::vector<double>* D0OpAngleMatching = 0;
	std::vector<double>* deltaRD0Matching = 0;

	double Variable_D0etaMatching = 0.;
	double Variable_MCD0etaMatching = 0.;
	double Variable_D0phiMatching = 0.;
	double Variable_MCD0phiMatching = 0.;
	double Variable_D0ptMatching = 0.;
	double Variable_MCD0ptMatching = 0.;
	double Variable_D0KtMatching = 0.;
	double Variable_D0SxyMatching = 0.;
	double Variable_D0OpAngleMatching = 0.;
	double Variable_deltaRD0Matching = 0.;

	string Variable_NameOfFiredTriggers; Variable_NameOfFiredTriggers.clear();
	
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
	//For D0 
	//-----------------------------------------
	t_D0analysis.Branch("D0Kpimass",&Variable_D0Kpimass, "Variable_D0Kpimass/D");
	t_D0analysis.Branch("D0Kpi_VtxProb",&Variable_D0Kpi_VtxProb, "Variable_D0Kpi_VtxProb/D");
	t_D0analysis.Branch("D0Kpis3D",&Variable_D0Kpis3D, "Variable_D0Kpis3D/D");
	t_D0analysis.Branch("D0KpiDispAngle",&Variable_D0KpiDispAngle, "Variable_D0KpiDispAngle/D");

	t_D0analysis.Branch("D0Kpipt",&Variable_D0Kpipt, "Variable_D0Kpipt/D");
	t_D0analysis.Branch("D0Kpieta",&Variable_D0Kpieta, "Variable_D0Kpieta/D");
	t_D0analysis.Branch("D0Kpiphi",&Variable_D0Kpiphi, "Variable_D0Kpiphi/D");

	t_D0analysis.Branch("TrkD0Kdxy",&Variable_TrkD0Kdxy, "Variable_TrkD0Kdxy/D");
	t_D0analysis.Branch("TrkD0Kdz",&Variable_TrkD0Kdz, "Variable_TrkD0Kdz/D");
	t_D0analysis.Branch("TrkD0Kchi2",&Variable_TrkD0Kchi2, "Variable_TrkD0Kchi2/D");
	t_D0analysis.Branch("TrkD0Kpt",&Variable_TrkD0Kpt, "Variable_TrkD0Kpt/D");
	t_D0analysis.Branch("TrkD0Keta",&Variable_TrkD0Keta, "Variable_TrkD0Keta/D");
	t_D0analysis.Branch("TrkD0kphi",&Variable_TrkD0kphi, "Variable_TrkD0kphi/D");
	t_D0analysis.Branch("TrkD0Knhits",&Variable_TrkD0Knhits, "Variable_TrkD0Knhits/D");

	t_D0analysis.Branch("TrkD0pidxy",&Variable_TrkD0pidxy, "Variable_TrkD0pidxy/D");
	t_D0analysis.Branch("TrkD0pidz",&Variable_TrkD0pidz, "Variable_TrkD0pidz/D");
	t_D0analysis.Branch("TrkD0pichi2",&Variable_TrkD0pichi2, "Variable_TrkD0pichi2/D");
	t_D0analysis.Branch("TrkD0pipt",&Variable_TrkD0pipt, "Variable_TrkD0pipt/D");
	t_D0analysis.Branch("TrkD0pieta",&Variable_TrkD0pieta, "Variable_TrkD0pieta/D");
	t_D0analysis.Branch("TrkD0piphi",&Variable_TrkD0piphi, "Variable_TrkD0piphi/D");
	t_D0analysis.Branch("TrkD0pinhits",&Variable_TrkD0pinhits, "Variable_TrkD0pinhits/D");

	t_D0analysis.Branch("D0KpisXY",&Variable_D0KpisXY, "Variable_D0KpisXY/D");
	t_D0analysis.Branch("D0Kpid3D",&Variable_D0Kpid3D, "Variable_D0Kpid3D/D");
	t_D0analysis.Branch("D0Kpie3D",&Variable_D0Kpie3D, "Variable_D0Kpie3D/D");
	t_D0analysis.Branch("D0_kT",&Variable_D0_kT, "Variable_D0_kT/D");

	//----------------------------------------
	//For D* MC
	//-----------------------------------------
	t_DsMC.Branch("dScandsKpi",&Variable_dScandsKpi, "Variable_dScandsKpi/D");
	t_DsMC.Branch("MCDseta",&Variable_MCDseta, "Variable_MCDseta/D");
	t_DsMC.Branch("MCDsphi",&Variable_MCDsphi, "Variable_MCDsphi/D");
	t_DsMC.Branch("MCDspt",&Variable_MCDspt, "Variable_MCDspt/D");
	t_DsMC.Branch("MCDsenergy",&Variable_MCDsenergy, "Variable_MCDsenergy/D");
	t_DsMC.Branch("MCDsp",&Variable_MCDsp, "Variable_MCDsp/D");
	t_DsMC.Branch("MCDset",&Variable_MCDset, "Variable_MCDset/D");
	t_DsMC.Branch("MCDsmass",&Variable_MCDsmass, "Variable_MCDsmass/D");
	t_DsMC.Branch("MCD0eta",&Variable_MCD0eta, "Variable_MCD0eta/D");
	t_DsMC.Branch("MCD0phi",&Variable_MCD0phi, "Variable_MCD0phi/D");
	t_DsMC.Branch("MCD0pt",&Variable_MCD0pt, "Variable_MCD0pt/D");
	t_DsMC.Branch("MCD0energy",&Variable_MCD0energy, "Variable_MCD0energy/D");
	t_DsMC.Branch("MCD0p",&Variable_MCD0p, "Variable_MCD0p/D");
	t_DsMC.Branch("MCD0et",&Variable_MCD0et, "Variable_MCD0et/D");
	t_DsMC.Branch("MCD0rapidity",&Variable_MCD0rapidity, "Variable_MCD0rapidity/D");
	t_DsMC.Branch("MCD0mass",&Variable_MCD0mass, "Variable_MCD0mass/D");
	t_DsMC.Branch("MCDsKphi",&Variable_MCDsKphi, "Variable_MCDsKphi/D");
	t_DsMC.Branch("MCDsKpt",&Variable_MCDsKpt, "Variable_MCDsKpt/D");
	t_DsMC.Branch("MCDsKenergy",&Variable_MCDsKenergy, "Variable_MCDsKenergy/D");
	t_DsMC.Branch("MCDsKp",&Variable_MCDsKp, "Variable_MCDsKp/D");
	t_DsMC.Branch("MCDsKet",&Variable_MCDsKet, "Variable_MCDsKet/D");
	t_DsMC.Branch("MCDsKrapidity",&Variable_MCDsKrapidity, "Variable_MCDsKrapidity/D");
	t_DsMC.Branch("MCDsKmass",&Variable_MCDsKmass, "Variable_MCDsKmass/D");
	t_DsMC.Branch("MCDsPieta",&Variable_MCDsPieta, "Variable_MCDsPieta/D");
	t_DsMC.Branch("MCDsPiphi",&Variable_MCDsPiphi, "Variable_MCDsPiphi/D");
	t_DsMC.Branch("MCDsPipt",&Variable_MCDsPipt, "Variable_MCDsPipt/D");
	t_DsMC.Branch("MCDsPienergy",&Variable_MCDsPienergy, "Variable_MCDsPienergy/D");
	t_DsMC.Branch("MCDsPip",&Variable_MCDsPip, "Variable_MCDsPip/D");
	t_DsMC.Branch("MCDsPiet",&Variable_MCDsPiet, "Variable_MCDsPiet/D");
	t_DsMC.Branch("MCDsPirapidity",&Variable_MCDsPirapidity, "Variable_MCDsPirapidity/D");
	t_DsMC.Branch("MCDsPimass",&Variable_MCDsPimass, "Variable_MCDsPimass/D");

	//----------------------------------------
	//For D0 MC
	//-----------------------------------------
	t_D0MC.Branch("MCpromptD0eta",&Variable_MCpromptD0eta, "Variable_MCpromptD0eta/D");
	t_D0MC.Branch("MCpromptD0phi",&Variable_MCpromptD0phi, "Variable_MCpromptD0phi/D");
	t_D0MC.Branch("MCpromptD0pt",&Variable_MCpromptD0pt, "Variable_MCpromptD0pt/D");
	t_D0MC.Branch("MCpromptD0energy",&Variable_MCpromptD0energy, "Variable_MCpromptD0energy/D");
	t_D0MC.Branch("MCpromptD0p",&Variable_MCpromptD0p, "Variable_MCpromptD0p/D");
	t_D0MC.Branch("MCpromptD0et",&Variable_MCpromptD0et, "Variable_MCpromptD0et/D");
	t_D0MC.Branch("MCpromptD0rapidity",&Variable_MCpromptD0rapidity, "Variable_MCpromptD0rapidity/D");
	t_D0MC.Branch("MCpromptD0mass",&Variable_MCpromptD0mass, "Variable_MCpromptD0mass/D");
	t_D0MC.Branch("MCpromptD0_Keta",&Variable_MCpromptD0_Keta, "Variable_MCpromptD0_Keta/D");
	t_D0MC.Branch("MCpromptD0_Kphi",&Variable_MCpromptD0_Kphi, "Variable_MCpromptD0_Kphi/D");
	t_D0MC.Branch("MCpromptD0_Kpt",&Variable_MCpromptD0_Kpt, "Variable_MCpromptD0_Kpt/D");
	t_D0MC.Branch("MCpromptD0_Kenergy",&Variable_MCpromptD0_Kenergy, "Variable_MCpromptD0_Kenergy/D");
	t_D0MC.Branch("MCpromptD0_Kp",&Variable_MCpromptD0_Kp, "Variable_MCpromptD0_Kp/D");
	t_D0MC.Branch("MCpromptD0_Ket",&Variable_MCpromptD0_Ket, "Variable_MCpromptD0_Ket/D");
	t_D0MC.Branch("MCpromptD0_Krapidity",&Variable_MCpromptD0_Krapidity, "Variable_MCpromptD0_Krapidity/D");
	t_D0MC.Branch("MCpromptD0_Kmass",&Variable_MCpromptD0_Kmass, "Variable_MCpromptD0_Kmass/D");
	t_D0MC.Branch("MCpromptD0_Pieta",&Variable_MCpromptD0_Pieta, "Variable_MCpromptD0_Pieta/D");
	t_D0MC.Branch("MCpromptD0_Piphi",&Variable_MCpromptD0_Piphi, "Variable_MCpromptD0_Piphi/D");
	t_D0MC.Branch("MCpromptD0_Pipt",&Variable_MCpromptD0_Pipt, "Variable_MCpromptD0_Pipt/D");
	t_D0MC.Branch("MCpromptD0_Pienergy",&Variable_MCpromptD0_Pienergy, "Variable_MCpromptD0_Pienergy/D");
	t_D0MC.Branch("MCpromptD0_Pip",&Variable_MCpromptD0_Pip, "Variable_MCpromptD0_Pip/D");
	t_D0MC.Branch("MCpromptD0_Piet",&Variable_MCpromptD0_Piet, "Variable_MCpromptD0_Piet/D");
	t_D0MC.Branch("MCpromptD0_Pirapidity",&Variable_MCpromptD0_Pirapidity, "Variable_MCpromptD0_Pirapidity/D");
	t_D0MC.Branch("MCpromptD0_Pimass",&Variable_MCpromptD0_Pimass, "Variable_MCpromptD0_Pimass/D");
	t_D0MC.Branch("MCpromptD0_DispAngle",&Variable_MCpromptD0_DispAngle, "Variable_MCpromptD0_DispAngle/D");

	//----------------------------------------
	//For D* Matching
	//-----------------------------------------
	t_DsMatching.Branch("DsetaMatching",&Variable_DsetaMatching, "Variable_DsetaMatching/D");
	t_DsMatching.Branch("MCDsetaMatching",&Variable_MCDsetaMatching, "Variable_MCDsetaMatching/D");
	t_DsMatching.Branch("DsphiMatching",&Variable_DsphiMatching, "Variable_DsphiMatching/D");
	t_DsMatching.Branch("MCDsphiMatching",&Variable_MCDsphiMatching, "Variable_MCDsphiMatching/D");
	t_DsMatching.Branch("DsptMatching",&Variable_DsptMatching, "Variable_DsptMatching/D");
	t_DsMatching.Branch("MCDsptMatching",&Variable_MCDsptMatching, "Variable_MCDsptMatching/D");
	t_DsMatching.Branch("D0fromDsmass",&Variable_D0fromDsmass, "Variable_D0fromDsmass/D");
	t_DsMatching.Branch("deltaRDsMatching",&Variable_deltaRDsMatching, "Variable_deltaRDsMatching/D");
	
	//----------------------------------------
	//For D0 Matching
	//-----------------------------------------
	t_D0Matching.Branch("D0etaMatching",&Variable_D0etaMatching, "Variable_D0etaMatching/D");
	t_D0Matching.Branch("MCD0etaMatching",&Variable_MCD0etaMatching, "Variable_MCD0etaMatching/D");
	t_D0Matching.Branch("D0phiMatching",&Variable_D0phiMatching, "Variable_D0phiMatching/D");
	t_D0Matching.Branch("MCD0phiMatching",&Variable_MCD0phiMatching, "Variable_MCD0phiMatching/D");
	t_D0Matching.Branch("D0ptMatching",&Variable_D0ptMatching, "Variable_D0ptMatching/D");
	t_D0Matching.Branch("MCD0ptMatching",&Variable_MCD0ptMatching, "Variable_MCD0ptMatching/D");
	t_D0Matching.Branch("D0KtMatching",&Variable_D0KtMatching, "Variable_D0KtMatching/D");
	t_D0Matching.Branch("D0SxyMatching",&Variable_D0SxyMatching, "Variable_D0SxyMatching/D");
	t_D0Matching.Branch("D0OpAngleMatching",&Variable_D0OpAngleMatching, "Variable_D0OpAngleMatching/D");
	t_D0Matching.Branch("deltaRD0Matching",&Variable_deltaRD0Matching, "Variable_deltaRD0Matching/D");


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

	TH1 *hRecEta_Rec1 = makeTH1("hRecEta_Rec1", 100,0,5, "hRecEta_Rec1 ; #eta ; Events ", kRed);
	TH1 *hRecEta_MC1 = makeTH1("hRecEta_MC1", 100,0,5, "hRecEta_MC1 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec1 = makeTH1("hRecPhi_Rec1", 100,0,5, "hRecPhi_Rec1 ; #phi ; Events ", kRed);
	TH1 *hRecPhi_MC1 = makeTH1("hRecPhi_MC1", 100,0,5, "hRecPhi_MC1 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec1 = makeTH1("hRecPt_Rec1", 100,0,20, "hRecPt_Rec1 ; pT ; Events ", kRed);
	TH1 *hRecPt_MC1 = makeTH1("hRecPt_MC1", 100,0,20, "hRecPt_MC1 ; pT ; Events ", kRed);

	TH1 *hRecEta_Rec2 = makeTH1("hRecEta_Rec2", 100,0,5, "hRecEta_Rec2 ; #eta ; Events ", kRed);
	TH1 *hRecEta_MC2 = makeTH1("hRecEta_MC2", 100,0,5, "hRecEta_MC2 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec2 = makeTH1("hRecPhi_Rec2", 100,0,5, "hRecPhi_Rec2 ; #phi ; Events ", kRed);
	TH1 *hRecPhi_MC2 = makeTH1("hRecPhi_MC2", 100,0,5, "hRecPhi_MC2 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec2 = makeTH1("hRecPt_Rec2", 100,0,20, "hRecPt_Rec2 ; pT ; Events ", kRed);
	TH1 *hRecPt_MC2 = makeTH1("hRecPt_MC2", 100,0,20, "hRecPt_MC2 ; pT ; Events ", kRed);

	TH1 *hRecEta_Rec3 = makeTH1("hRecEta_Rec3", 100,0,5, "hRecEta_Rec3 ; #eta ; Events ", kRed);
	TH1 *hRecEta_MC3 = makeTH1("hRecEta_MC3", 100,0,5, "hRecEta_MC3 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec3 = makeTH1("hRecPhi_Rec3", 100,0,5, "hRecPhi_Rec3 ; #phi ; Events ", kRed);
	TH1 *hRecPhi_MC3 = makeTH1("hRecPhi_MC3", 100,0,5, "hRecPhi_MC3 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec3 = makeTH1("hRecPt_Rec3", 100,0,20, "hRecPt_Rec3 ; pT ; Events ", kRed);
	TH1 *hRecPt_MC3 = makeTH1("hRecPt_MC3", 100,0,20, "hRecPt_MC3 ; pT ; Events ", kRed);

	TH1 *hRecEta_Rec4 = makeTH1("hRecEta_Rec4", 100,0,5, "hRecEta_Rec4 ; #eta ; Events ", kRed);
	TH1 *hRecEta_MC4 = makeTH1("hRecEta_MC4", 100,0,5, "hRecEta_MC4 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec4 = makeTH1("hRecPhi_Rec4", 100,0,5, "hRecPhi_Rec4 ; #phi ; Events ", kRed);
	TH1 *hRecPhi_MC4 = makeTH1("hRecPhi_MC4", 100,0,5, "hRecPhi_MC4 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec4 = makeTH1("hRecPt_Rec4", 100,0,20, "hRecPt_Rec4 ; pT ; Events ", kRed);
	TH1 *hRecPt_MC4 = makeTH1("hRecPt_MC4", 100,0,20, "hRecPt_MC4 ; pT ; Events ", kRed);

	TH1 *hRecEta_Rec5 = makeTH1("hRecEta_Rec5", 100,0,5, "hRecEta_Rec5 ; #eta ; Events ", kRed);
	TH1 *hRecEta_MC5 = makeTH1("hRecEta_MC5", 100,0,5, "hRecEta_MC5 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec5 = makeTH1("hRecPhi_Rec5", 100,0,5, "hRecPhi_Rec5 ; #phi ; Events ", kRed);
	TH1 *hRecPhi_MC5 = makeTH1("hRecPhi_MC5", 100,0,5, "hRecPhi_MC5 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec5 = makeTH1("hRecPt_Rec5", 100,0,20, "hRecPt_Rec5 ; pT ; Events ", kRed);
	TH1 *hRecPt_MC5 = makeTH1("hRecPt_MC5", 100,0,20, "hRecPt_MC5 ; pT ; Events ", kRed);

	TH1 *hRecEta_Rec6 = makeTH1("hRecEta_Rec6", 100,0,5, "hRecEta_Rec6 ; #eta ; Events ", kRed);
	TH1 *hRecEta_MC6 = makeTH1("hRecEta_MC6", 100,0,5, "hRecEta_MC6 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec6 = makeTH1("hRecPhi_Rec6", 100,0,5, "hRecPhi_Rec6 ; #phi ; Events ", kRed);
	TH1 *hRecPhi_MC6 = makeTH1("hRecPhi_MC6", 100,0,5, "hRecPhi_MC6 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec6 = makeTH1("hRecPt_Rec6", 100,0,20, "hRecPt_Rec6 ; pT ; Events ", kRed);
	TH1 *hRecPt_MC6 = makeTH1("hRecPt_MC6", 100,0,20, "hRecPt_MC6 ; pT ; Events ", kRed);

	//for D0
	TH1 *hRecD0Eta_Rec1 = makeTH1("hRecD0Eta_Rec1", 100,0,5, "hRecD0Eta_Rec1 ; #eta ; Events ", kRed);
	TH1 *hRecD0Eta_MC1 = makeTH1("hRecD0Eta_MC1", 100,0,5, "hRecD0Eta_MC1 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec1 = makeTH1("hRecD0Phi_Rec1", 100,0,5, "hRecD0Phi_Rec1 ; #phi ; Events ", kRed);
	TH1 *hRecD0Phi_MC1 = makeTH1("hRecD0Phi_MC1", 100,0,5, "hRecD0Phi_MC1 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec1 = makeTH1("hRecD0Pt_Rec1", 100,0,20, "hRecD0Pt_Rec1 ; pT ; Events ", kRed);
	TH1 *hRecD0Pt_MC1 = makeTH1("hRecD0Pt_MC1", 100,0,20, "hRecD0Pt_MC1 ; pT ; Events ", kRed);

	TH1 *hRecD0Eta_Rec2 = makeTH1("hRecD0Eta_Rec2", 100,0,5, "hRecD0Eta_Rec2 ; #eta ; Events ", kRed);
	TH1 *hRecD0Eta_MC2 = makeTH1("hRecD0Eta_MC2", 100,0,5, "hRecD0Eta_MC2 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec2 = makeTH1("hRecD0Phi_Rec2", 100,0,5, "hRecD0Phi_Rec2 ; #phi ; Events ", kRed);
	TH1 *hRecD0Phi_MC2 = makeTH1("hRecD0Phi_MC2", 100,0,5, "hRecD0Phi_MC2 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec2 = makeTH1("hRecD0Pt_Rec2", 100,0,20, "hRecD0Pt_Rec2 ; pT ; Events ", kRed);
	TH1 *hRecD0Pt_MC2 = makeTH1("hRecD0Pt_MC2", 100,0,20, "hRecD0Pt_MC2 ; pT ; Events ", kRed);
	
	TH1 *hRecD0Eta_Rec3 = makeTH1("hRecD0Eta_Rec3", 100,0,5, "hRecD0Eta_Rec3 ; #eta ; Events ", kRed);
	TH1 *hRecD0Eta_MC3 = makeTH1("hRecD0Eta_MC3", 100,0,5, "hRecD0Eta_MC3 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec3 = makeTH1("hRecD0Phi_Rec3", 100,0,5, "hRecD0Phi_Rec3 ; #phi ; Events ", kRed);
	TH1 *hRecD0Phi_MC3 = makeTH1("hRecD0Phi_MC3", 100,0,5, "hRecD0Phi_MC3 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec3 = makeTH1("hRecD0Pt_Rec3", 100,0,20, "hRecD0Pt_Rec3 ; pT ; Events ", kRed);
	TH1 *hRecD0Pt_MC3 = makeTH1("hRecD0Pt_MC3", 100,0,20, "hRecD0Pt_MC3 ; pT ; Events ", kRed);

	TH1 *hRecD0Eta_Rec4 = makeTH1("hRecD0Eta_Rec4", 100,0,5, "hRecD0Eta_Rec4 ; #eta ; Events ", kRed);
	TH1 *hRecD0Eta_MC4 = makeTH1("hRecD0Eta_MC4", 100,0,5, "hRecD0Eta_MC4 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec4 = makeTH1("hRecD0Phi_Rec4", 100,0,5, "hRecD0Phi_Rec4 ; #phi ; Events ", kRed);
	TH1 *hRecD0Phi_MC4 = makeTH1("hRecD0Phi_MC4", 100,0,5, "hRecD0Phi_MC4 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec4 = makeTH1("hRecD0Pt_Rec4", 100,0,20, "hRecD0Pt_Rec4 ; pT ; Events ", kRed);
	TH1 *hRecD0Pt_MC4 = makeTH1("hRecD0Pt_MC4", 100,0,20, "hRecD0Pt_MC4 ; pT ; Events ", kRed);

	TH1 *hRecD0Eta_Rec5 = makeTH1("hRecD0Eta_Rec5", 100,0,5, "hRecD0Eta_Rec5 ; #eta ; Events ", kRed);
	TH1 *hRecD0Eta_MC5 = makeTH1("hRecD0Eta_MC5", 100,0,5, "hRecD0Eta_MC5 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec5 = makeTH1("hRecD0Phi_Rec5", 100,0,5, "hRecD0Phi_Rec5 ; #phi ; Events ", kRed);
	TH1 *hRecD0Phi_MC5 = makeTH1("hRecD0Phi_MC5", 100,0,5, "hRecD0Phi_MC5 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec5 = makeTH1("hRecD0Pt_Rec5", 100,0,20, "hRecD0Pt_Rec5 ; pT ; Events ", kRed);
	TH1 *hRecD0Pt_MC5 = makeTH1("hRecD0Pt_MC5", 100,0,20, "hRecD0Pt_MC5 ; pT ; Events ", kRed);

	TH1 *hRecD0Eta_Rec6 = makeTH1("hRecD0Eta_Rec6", 100,0,5, "hRecD0Eta_Rec6 ; #eta ; Events ", kRed);
	TH1 *hRecD0Eta_MC6 = makeTH1("hRecD0Eta_MC6", 100,0,5, "hRecD0Eta_MC6 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec6 = makeTH1("hRecD0Phi_Rec6", 100,0,5, "hRecD0Phi_Rec6 ; #phi ; Events ", kRed);
	TH1 *hRecD0Phi_MC6 = makeTH1("hRecD0Phi_MC6", 100,0,5, "hRecD0Phi_MC6 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec6 = makeTH1("hRecD0Pt_Rec6", 100,0,20, "hRecD0Pt_Rec6 ; pT ; Events ", kRed);
	TH1 *hRecD0Pt_MC6 = makeTH1("hRecD0Pt_MC6", 100,0,20, "hRecD0Pt_MC6 ; pT ; Events ", kRed);

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
	//FOR D0
	//-----------------------------------------
	TBranch *b_D0Kpimass ; t1->SetBranchAddress("D0Kpimass",&D0Kpimass,&b_D0Kpimass);
	TBranch *b_D0Kpi_VtxProb ; t1->SetBranchAddress("D0Kpi_VtxProb",&D0Kpi_VtxProb,&b_D0Kpi_VtxProb);
	TBranch *b_D0Kpis3D ; t1->SetBranchAddress("D0Kpis3D",&D0Kpis3D,&b_D0Kpis3D);
	TBranch *b_D0KpiDispAngle ; t1->SetBranchAddress("D0KpiDispAngle",&D0KpiDispAngle,&b_D0KpiDispAngle);

	TBranch *b_D0Kpipt ; t1->SetBranchAddress("D0Kpipt",&D0Kpipt,&b_D0Kpipt);
	TBranch *b_D0Kpieta ; t1->SetBranchAddress("D0Kpieta",&D0Kpieta,&b_D0Kpieta);
	TBranch *b_D0Kpiphie ; t1->SetBranchAddress("D0Kpiphi",&D0Kpiphi,&b_D0Kpiphie);

	TBranch *b_TrkD0Kdxy ; t1->SetBranchAddress("TrkD0Kdxy",&TrkD0Kdxy,&b_TrkD0Kdxy);
	TBranch *b_TrkD0Kdz ; t1->SetBranchAddress("TrkD0Kdz",&TrkD0Kdz,&b_TrkD0Kdz);
	TBranch *b_TrkD0Kchi2 ; t1->SetBranchAddress("TrkD0Kchi2",&TrkD0Kchi2,&b_TrkD0Kchi2);
	TBranch *b_TrkD0Kpt ; t1->SetBranchAddress("TrkD0Kpt",&TrkD0Kpt,&b_TrkD0Kpt);
	TBranch *b_TrkD0Keta ; t1->SetBranchAddress("TrkD0Keta",&TrkD0Keta,&b_TrkD0Keta);
	TBranch *b_TrkD0kphi ; t1->SetBranchAddress("TrkD0kphi",&TrkD0kphi,&b_TrkD0kphi);
	TBranch *b_TrkD0Knhits ; t1->SetBranchAddress("TrkD0Knhits",&TrkD0Knhits,&b_TrkD0Knhits);

	TBranch *b_TrkD0pidxy ; t1->SetBranchAddress("TrkD0pidxy",&TrkD0pidxy,&b_TrkD0pidxy);
	TBranch *b_TrkD0pidz ; t1->SetBranchAddress("TrkD0pidz",&TrkD0pidz,&b_TrkD0pidz);
	TBranch *b_TrkD0pichi2 ; t1->SetBranchAddress("TrkD0pichi2",&TrkD0pichi2,&b_TrkD0pichi2);
	TBranch *b_TrkD0pipt ; t1->SetBranchAddress("TrkD0pipt",&TrkD0pipt,&b_TrkD0pipt);
	TBranch *b_TrkD0pieta ; t1->SetBranchAddress("TrkD0pieta",&TrkD0pieta,&b_TrkD0pieta);
	TBranch *b_TrkD0piphi ; t1->SetBranchAddress("TrkD0piphi",&TrkD0piphi,&b_TrkD0piphi);
	TBranch *b_TrkD0pinhits ; t1->SetBranchAddress("TrkD0pinhits",&TrkD0pinhits,&b_TrkD0pinhits);

	TBranch *b_D0KpisXY ; t1->SetBranchAddress("D0KpisXY",&D0KpisXY,&b_D0KpisXY);
	TBranch *b_D0Kpid3D ; t1->SetBranchAddress("D0Kpid3D",&D0Kpid3D,&b_D0Kpid3D);
	TBranch *b_D0Kpie3D ; t1->SetBranchAddress("D0Kpie3D",&D0Kpie3D,&b_D0Kpie3D);
	TBranch *b_D0_kT ; t1->SetBranchAddress("D0KpikT",&D0_kT,&b_D0_kT);

	//----------------------------------------
	//FOR D* MC
	//-----------------------------------------
	TBranch *b_dScandsKpi ; t1->SetBranchAddress("dScandsKpi",&dScandsKpi,&b_dScandsKpi);
	TBranch *b_MCDseta ; t1->SetBranchAddress("MCDseta",&MCDseta,&b_MCDseta);
	TBranch *b_MCDsphi ; t1->SetBranchAddress("MCDsphi",&MCDsphi,&b_MCDsphi);
	TBranch *b_MCDspt ; t1->SetBranchAddress("MCDspt",&MCDspt,&b_MCDspt);
	TBranch *b_MCDsenergy ; t1->SetBranchAddress("MCDsenergy",&MCDsenergy,&b_MCDsenergy);
	TBranch *b_MCDsp ; t1->SetBranchAddress("MCDsp",&MCDsp,&b_MCDsp);
	TBranch *b_MCDset ; t1->SetBranchAddress("MCDset",&MCDset,&b_MCDset);
	TBranch *b_MCDsmass ; t1->SetBranchAddress("MCDsmass",&MCDsmass,&b_MCDsmass);
	TBranch *b_MCD0eta ; t1->SetBranchAddress("MCD0eta",&MCD0eta,&b_MCD0eta);
	TBranch *b_MCD0phi ; t1->SetBranchAddress("MCD0phi",&MCD0phi,&b_MCD0phi);
	TBranch *b_MCD0pt ; t1->SetBranchAddress("MCD0pt",&MCD0pt,&b_MCD0pt);
	TBranch *b_MCD0energy ; t1->SetBranchAddress("MCD0energy",&MCD0energy,&b_MCD0energy);
	TBranch *b_MCD0p ; t1->SetBranchAddress("MCD0p",&MCD0p,&b_MCD0p);
	TBranch *b_MCD0et ; t1->SetBranchAddress("MCD0et",&MCD0et,&b_MCD0et);
	TBranch *b_MCD0rapidity ; t1->SetBranchAddress("MCD0rapidity",&MCD0rapidity,&b_MCD0rapidity);
	TBranch *b_MCD0mass ; t1->SetBranchAddress("MCD0mass",&MCD0mass,&b_MCD0mass);
	TBranch *b_MCDsKphi ; t1->SetBranchAddress("MCDsKphi",&MCDsKphi,&b_MCDsKphi);
	TBranch *b_MCDsKpt ; t1->SetBranchAddress("MCDsKpt",&MCDsKpt,&b_MCDsKpt);
	TBranch *b_MCDsKenergy ; t1->SetBranchAddress("MCDsKenergy",&MCDsKenergy,&b_MCDsKenergy);
	TBranch *b_MCDsKp ; t1->SetBranchAddress("MCDsKp",&MCDsKp,&b_MCDsKp);
	TBranch *b_MCDsKet ; t1->SetBranchAddress("MCDsKet",&MCDsKet,&b_MCDsKet);
	TBranch *b_MCDsKrapidity ; t1->SetBranchAddress("MCDsKrapidity",&MCDsKrapidity,&b_MCDsKrapidity);
	TBranch *b_MCDsKmass ; t1->SetBranchAddress("MCDsKmass",&MCDsKmass,&b_MCDsKmass);
	TBranch *b_MCDsPieta ; t1->SetBranchAddress("MCDsPieta",&MCDsPieta,&b_MCDsPieta);
	TBranch *b_MCDsPiphi ; t1->SetBranchAddress("MCDsPiphi",&MCDsPiphi,&b_MCDsPiphi);
	TBranch *b_MCDsPipt ; t1->SetBranchAddress("MCDsPipt",&MCDsPipt,&b_MCDsPipt);
	TBranch *b_MCDsPienergy ; t1->SetBranchAddress("MCDsPienergy",&MCDsPienergy,&b_MCDsPienergy);
	TBranch *b_MCDsPip ; t1->SetBranchAddress("MCDsPip",&MCDsPip,&b_MCDsPip);
	TBranch *b_MCDsPiet ; t1->SetBranchAddress("MCDsPiet",&MCDsPiet,&b_MCDsPiet);
	TBranch *b_MCDsPirapidity ; t1->SetBranchAddress("MCDsPirapidity",&MCDsPirapidity,&b_MCDsPirapidity);
	TBranch *b_MCDsPimass ; t1->SetBranchAddress("MCDsPimass",&MCDsPimass,&b_MCDsPimass);

	//----------------------------------------
	//FOR D0 MC
	//-----------------------------------------
	TBranch *b_MCpromptD0eta ; t1->SetBranchAddress("MCpromptD0eta",&MCpromptD0eta,&b_MCpromptD0eta);
	TBranch *b_MCpromptD0phi ; t1->SetBranchAddress("MCpromptD0phi",&MCpromptD0phi,&b_MCpromptD0phi);
	TBranch *b_MCpromptD0pt ; t1->SetBranchAddress("MCpromptD0pt",&MCpromptD0pt,&b_MCpromptD0pt);
	TBranch *b_MCpromptD0energy ; t1->SetBranchAddress("MCpromptD0energy",&MCpromptD0energy,&b_MCpromptD0energy);
	TBranch *b_MCpromptD0p ; t1->SetBranchAddress("MCpromptD0p",&MCpromptD0p,&b_MCpromptD0p);
	TBranch *b_MCpromptD0et ; t1->SetBranchAddress("MCpromptD0et",&MCpromptD0et,&b_MCpromptD0et);
	TBranch *b_MCpromptD0rapidity ; t1->SetBranchAddress("MCpromptD0rapidity",&MCpromptD0rapidity,&b_MCpromptD0rapidity);
	TBranch *b_MCpromptD0mass ; t1->SetBranchAddress("MCpromptD0mass",&MCpromptD0mass,&b_MCpromptD0mass);
	TBranch *b_MCpromptD0_Keta ; t1->SetBranchAddress("MCpromptD0_Keta",&MCpromptD0_Keta,&b_MCpromptD0_Keta);
	TBranch *b_MCpromptD0_Kphi ; t1->SetBranchAddress("MCpromptD0_Kphi",&MCpromptD0_Kphi,&b_MCpromptD0_Kphi);
	TBranch *b_MCpromptD0_Kpt ; t1->SetBranchAddress("MCpromptD0_Kpt",&MCpromptD0_Kpt,&b_MCpromptD0_Kpt);
	TBranch *b_MCpromptD0_Kenergy ; t1->SetBranchAddress("MCpromptD0_Kenergy",&MCpromptD0_Kenergy,&b_MCpromptD0_Kenergy);
	TBranch *b_MCpromptD0_Kp ; t1->SetBranchAddress("MCpromptD0_Kp",&MCpromptD0_Kp,&b_MCpromptD0_Kp);
	TBranch *b_MCpromptD0_Ket ; t1->SetBranchAddress("MCpromptD0_Ket",&MCpromptD0_Ket,&b_MCpromptD0_Ket);
	TBranch *b_MCpromptD0_Krapidity ; t1->SetBranchAddress("MCpromptD0_Krapidity",&MCpromptD0_Krapidity,&b_MCpromptD0_Krapidity);
	TBranch *b_MCpromptD0_Kmass ; t1->SetBranchAddress("MCpromptD0_Kmass",&MCpromptD0_Kmass,&b_MCpromptD0_Kmass);
	TBranch *b_MCpromptD0_Pieta ; t1->SetBranchAddress("MCpromptD0_Pieta",&MCpromptD0_Pieta,&b_MCpromptD0_Pieta);
	TBranch *b_MCpromptD0_Piphi ; t1->SetBranchAddress("MCpromptD0_Piphi",&MCpromptD0_Piphi,&b_MCpromptD0_Piphi);
	TBranch *b_MCpromptD0_Pipt ; t1->SetBranchAddress("MCpromptD0_Pipt",&MCpromptD0_Pipt,&b_MCpromptD0_Pipt);
	TBranch *b_MCpromptD0_Pienergy ; t1->SetBranchAddress("MCpromptD0_Pienergy",&MCpromptD0_Pienergy,&b_MCpromptD0_Pienergy);
	TBranch *b_MCpromptD0_Pip ; t1->SetBranchAddress("MCpromptD0_Pip",&MCpromptD0_Pip,&b_MCpromptD0_Pip);
	TBranch *b_MCpromptD0_Piet ; t1->SetBranchAddress("MCpromptD0_Piet",&MCpromptD0_Piet,&b_MCpromptD0_Piet);
	TBranch *b_MCpromptD0_Pirapidity ; t1->SetBranchAddress("MCpromptD0_Pirapidity",&MCpromptD0_Pirapidity,&b_MCpromptD0_Pirapidity);
	TBranch *b_MCpromptD0_Pimass ; t1->SetBranchAddress("MCpromptD0_Pimass",&MCpromptD0_Pimass,&b_MCpromptD0_Pimass);
	TBranch *b_MCpromptD0_DispAngle ; t1->SetBranchAddress("MCpromptD0_DispAngle",&MCpromptD0_DispAngle,&b_MCpromptD0_DispAngle);

	//----------------------------------------
	//FOR D* Matching
	//-----------------------------------------
	TBranch *b_DsetaMatching ; t1->SetBranchAddress("DsetaMatching",&DsetaMatching,&b_DsetaMatching);
	TBranch *b_MCDsetaMatching ; t1->SetBranchAddress("MCDsetaMatching",&MCDsetaMatching,&b_MCDsetaMatching);
	TBranch *b_DsphiMatching ; t1->SetBranchAddress("DsphiMatching",&DsphiMatching,&b_DsphiMatching);
	TBranch *b_MCDsphiMatching ; t1->SetBranchAddress("MCDsphiMatching",&MCDsphiMatching,&b_MCDsphiMatching);
	TBranch *b_DsptMatching ; t1->SetBranchAddress("DsptMatching",&DsptMatching,&b_DsptMatching);
	TBranch *b_MCDsptMatching ; t1->SetBranchAddress("MCDsptMatching",&MCDsptMatching,&b_MCDsptMatching);
	TBranch *b_D0fromDsmass ; t1->SetBranchAddress("D0fromDsmass",&D0fromDsmass,&b_D0fromDsmass);
	
	TBranch *b_deltaRDsMatching ; t1->SetBranchAddress("deltaRDsMatching",&deltaRDsMatching,&b_deltaRDsMatching);

	//----------------------------------------
	//FOR D0 Matching
	//-----------------------------------------
	TBranch *b_D0etaMatching ; t1->SetBranchAddress("D0etaMatching",&D0etaMatching,&b_D0etaMatching);
	TBranch *b_MCD0etaMatching ; t1->SetBranchAddress("MCD0etaMatching",&MCD0etaMatching,&b_MCD0etaMatching);
	TBranch *b_D0phiMatching ; t1->SetBranchAddress("D0phiMatching",&D0phiMatching,&b_D0phiMatching);
	TBranch *b_MCD0phiMatching ; t1->SetBranchAddress("MCD0phiMatching",&MCD0phiMatching,&b_MCD0phiMatching);
	TBranch *b_D0ptMatching ; t1->SetBranchAddress("D0ptMatching",&D0ptMatching,&b_D0ptMatching);
	TBranch *b_MCD0ptMatching ; t1->SetBranchAddress("MCD0ptMatching",&MCD0ptMatching,&b_MCD0ptMatching);
	TBranch *b_D0KtMatching ; t1->SetBranchAddress("D0KtMatching",&D0KtMatching,&b_D0KtMatching);
	TBranch *b_D0SxyMatching ; t1->SetBranchAddress("D0SxyMatching",&D0SxyMatching,&b_D0SxyMatching);
	TBranch *b_D0OpAngleMatching ; t1->SetBranchAddress("D0OpAngleMatching",&D0OpAngleMatching,&b_D0OpAngleMatching);
	TBranch *b_deltaRD0Matching ; t1->SetBranchAddress("deltaRD0Matching",&deltaRD0Matching,&b_deltaRD0Matching);

	
	
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
		
		//----------------------------------------
		//To select the trigger
		//----------------------------------------
		b_NameOfFiredTriggers->GetEntry(ientry);
		bool goodtrigger = false;
		for(unsigned int tescont2=0; tescont2 < NameOfFiredTriggers->size(); tescont2++)
		{  string teststring = NameOfFiredTriggers->at(tescont2);
			if (teststring.find(Form("%s", TriggerPath.c_str() )) == 0) {goodtrigger = true;}
		}
		if(!goodtrigger) continue;
		//for(unsigned int tescont=0; tescont < NameOfFiredTriggers->size(); tescont++)
		//{  string teststring = NameOfFiredTriggers->at(tescont);
		//	std::cout << "NameOfFiredTriggers ["<< tescont <<"]: " << teststring  << std::endl;}
		
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
		//For D0 
		//-----------------------------------------
		b_D0Kpimass->GetEntry(ientry);
		b_D0Kpi_VtxProb->GetEntry(ientry);
		b_D0Kpis3D->GetEntry(ientry);
		b_D0KpiDispAngle->GetEntry(ientry);

		b_D0Kpipt->GetEntry(ientry);
		b_D0Kpieta->GetEntry(ientry);
		b_D0Kpiphie->GetEntry(ientry);
	
		b_TrkD0Kdxy->GetEntry(ientry);
		b_TrkD0Kdz->GetEntry(ientry);
		b_TrkD0Kchi2->GetEntry(ientry);
		b_TrkD0Kpt->GetEntry(ientry);
		b_TrkD0Keta->GetEntry(ientry);
		b_TrkD0kphi->GetEntry(ientry);
		b_TrkD0Knhits->GetEntry(ientry);

		b_TrkD0pidxy->GetEntry(ientry);
		b_TrkD0pidz->GetEntry(ientry);
		b_TrkD0pichi2->GetEntry(ientry);
		b_TrkD0pipt->GetEntry(ientry);
		b_TrkD0pieta->GetEntry(ientry);
		b_TrkD0piphi->GetEntry(ientry);
		b_TrkD0pinhits->GetEntry(ientry);

		b_D0KpisXY->GetEntry(ientry);
		b_D0Kpid3D->GetEntry(ientry);
		b_D0Kpie3D->GetEntry(ientry);
		b_D0_kT->GetEntry(ientry);

		//----------------------------------------
		//For D* MC
		//-----------------------------------------
		b_dScandsKpi->GetEntry(ientry);
		b_MCDseta->GetEntry(ientry);
		b_MCDsphi->GetEntry(ientry);
		b_MCDspt->GetEntry(ientry);
		b_MCDsenergy->GetEntry(ientry);
		b_MCDsp->GetEntry(ientry);
		b_MCDset->GetEntry(ientry);
		b_MCDsmass->GetEntry(ientry);
		b_MCD0eta->GetEntry(ientry);
		b_MCD0phi->GetEntry(ientry);
		b_MCD0pt->GetEntry(ientry);
		b_MCD0energy->GetEntry(ientry);
		b_MCD0p->GetEntry(ientry);
		b_MCD0et->GetEntry(ientry);
		b_MCD0rapidity->GetEntry(ientry);
		b_MCD0mass->GetEntry(ientry);
		b_MCD0mass->GetEntry(ientry);
		b_MCDsKphi->GetEntry(ientry);
		b_MCDsKpt->GetEntry(ientry);
		b_MCDsKenergy->GetEntry(ientry);
		b_MCDsKp->GetEntry(ientry);
		b_MCDsKet->GetEntry(ientry);
		b_MCDsKrapidity->GetEntry(ientry);
		b_MCDsKmass->GetEntry(ientry);
		b_MCDsPieta->GetEntry(ientry);
		b_MCDsPiphi->GetEntry(ientry);
		b_MCDsPipt->GetEntry(ientry);
		b_MCDsPienergy->GetEntry(ientry);
		b_MCDsPip->GetEntry(ientry);
		b_MCDsPiet->GetEntry(ientry);
		b_MCDsPirapidity->GetEntry(ientry);
		b_MCDsPimass->GetEntry(ientry);
	
		//----------------------------------------
		//For D0 MC
		//-----------------------------------------
		b_MCpromptD0eta->GetEntry(ientry);
		b_MCpromptD0phi->GetEntry(ientry);
		b_MCpromptD0pt->GetEntry(ientry);
		b_MCpromptD0energy->GetEntry(ientry);
		b_MCpromptD0p->GetEntry(ientry);
		b_MCpromptD0et->GetEntry(ientry);
		b_MCpromptD0rapidity->GetEntry(ientry);
		b_MCpromptD0mass->GetEntry(ientry);
		b_MCpromptD0_Keta->GetEntry(ientry);
		b_MCpromptD0_Kphi->GetEntry(ientry);
		b_MCpromptD0_Kpt->GetEntry(ientry);
		b_MCpromptD0_Kenergy->GetEntry(ientry);
		b_MCpromptD0_Kp->GetEntry(ientry);
		b_MCpromptD0_Ket->GetEntry(ientry);
		b_MCpromptD0_Krapidity->GetEntry(ientry);
		b_MCpromptD0_Kmass->GetEntry(ientry);
		b_MCpromptD0_Pieta->GetEntry(ientry);
		b_MCpromptD0_Piphi->GetEntry(ientry);
		b_MCpromptD0_Pipt->GetEntry(ientry);
		b_MCpromptD0_Pienergy->GetEntry(ientry);
		b_MCpromptD0_Pip->GetEntry(ientry);
		b_MCpromptD0_Piet->GetEntry(ientry);
		b_MCpromptD0_Pirapidity->GetEntry(ientry);
		b_MCpromptD0_Pimass->GetEntry(ientry);
		b_MCpromptD0_DispAngle->GetEntry(ientry);
		
		//----------------------------------------
		//For D* Matching 
		//-----------------------------------------
		b_DsetaMatching->GetEntry(ientry);
		b_MCDsetaMatching->GetEntry(ientry);
		b_DsphiMatching->GetEntry(ientry);
		b_MCDsphiMatching->GetEntry(ientry);
		b_DsptMatching->GetEntry(ientry);
		b_MCDsptMatching->GetEntry(ientry); 
		b_D0fromDsmass->GetEntry(ientry);
		b_deltaRDsMatching->GetEntry(ientry);

		//----------------------------------------
		//For D0 Matching
		//-----------------------------------------
		b_D0etaMatching->GetEntry(ientry);
		b_MCD0etaMatching->GetEntry(ientry);
		b_D0phiMatching->GetEntry(ientry);
		b_MCD0phiMatching->GetEntry(ientry);
		b_D0ptMatching->GetEntry(ientry);
		b_MCD0ptMatching->GetEntry(ientry);
		b_D0KtMatching->GetEntry(ientry);
		b_D0SxyMatching->GetEntry(ientry);
		b_D0OpAngleMatching->GetEntry(ientry);
		b_deltaRD0Matching->GetEntry(ientry);

		if (debug)cout << "debug 10 --------------------"<< endl;
	
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
			if (debug)cout << "debug D0(Direct) 12 --------------------"<< endl;
			if( D0Kpis3D->at(j) < SigCutD0) continue;
			if (debug)cout << "debug D0(Direct) 13 --------------------"<< endl;
			if( D0KpiDispAngle->at(j) < 0.99 ) continue;
			if (debug)cout << "debug D0(Direct) 14 --------------------"<< endl;
			//if ( D0Kpipt->at(i) < 3. ) continue;	
		
			D0KpimassHisto->Fill(D0Kpimass->at(j));
			D0Kpis3DHisto->Fill(D0Kpis3D->at(j));
			D0Kpi_VtxProbHisto->Fill(D0Kpi_VtxProb->at(j));

			Variable_D0Kpimass = D0Kpimass->at(j);
			Variable_D0Kpi_VtxProb = D0Kpi_VtxProb->at(j);
			Variable_D0Kpis3D = D0Kpis3D->at(j);
			Variable_D0KpiDispAngle = D0KpiDispAngle->at(j);

			Variable_D0Kpipt = D0Kpipt->at(j);
			Variable_D0Kpieta = D0Kpieta->at(j);
			Variable_D0Kpis3D = D0Kpis3D->at(j);
			Variable_D0Kpiphi = D0Kpiphi->at(j);

			Variable_TrkD0Kdxy = TrkD0Kdxy->at(j);
			Variable_TrkD0Kdz = TrkD0Kdz->at(j);
			Variable_TrkD0Kchi2 = TrkD0Kchi2->at(j);
			Variable_TrkD0Kpt = TrkD0Kpt->at(j);
			Variable_TrkD0Keta = TrkD0Keta->at(j);
			Variable_TrkD0kphi = TrkD0kphi->at(j);
			Variable_TrkD0Knhits = TrkD0Knhits->at(j);

			Variable_TrkD0pidxy = TrkD0pidxy->at(j);
			Variable_TrkD0pidz = TrkD0pidz->at(j);
			Variable_TrkD0pichi2 = TrkD0pichi2->at(j);
			Variable_TrkD0pieta = TrkD0pieta->at(j);
			Variable_TrkD0piphi = TrkD0piphi->at(j);
			Variable_TrkD0pinhits = TrkD0pinhits->at(j);

			Variable_D0KpisXY = D0KpisXY->at(j);
			Variable_D0Kpid3D = D0Kpid3D->at(j);
			Variable_D0Kpie3D = D0Kpie3D->at(j);
			Variable_D0_kT = D0_kT->at(j);
			
		
			//vectorInvariantMass_D0.push_back(D0mass->at(i));
			if (debug)cout << "debug D0(Direct) 12 --------------------"<< endl;
			t_D0analysis.Fill();	
  		}//For D0(Direct)

		//----------------------------------------
		//For D* MC
		//-----------------------------------------	
		for(unsigned int i=0; i < MCDsmass->size(); i++)
		{  
			//Variable_dScandsKpi = dScandsKpi->at(i);
			Variable_MCDseta = MCDseta->at(i);
			Variable_MCDsphi = MCDsphi->at(i);
			Variable_MCDspt = MCDspt->at(i);
			Variable_MCDsenergy = MCDsenergy->at(i);
			Variable_MCDsp = MCDsp->at(i);
			Variable_MCDset = MCDset->at(i);
			Variable_MCDsmass = MCDsmass->at(i);
			Variable_MCD0eta = MCD0eta->at(i);
			Variable_MCD0phi = MCD0phi->at(i);
			Variable_MCD0pt = MCD0pt->at(i);
			Variable_MCD0energy = MCD0energy->at(i);
			Variable_MCD0p = MCD0p->at(i);
			Variable_MCD0et = MCD0et->at(i);
			Variable_MCD0rapidity = MCD0rapidity->at(i);
			Variable_MCD0mass = MCD0mass->at(i);
			Variable_MCDsKphi = MCDsKphi->at(i);
			Variable_MCDsKpt = MCDsKpt->at(i);
			Variable_MCDsKenergy = MCDsKenergy->at(i);
			Variable_MCDsKp = MCDsKp->at(i);
			Variable_MCDsKet = MCDsKet->at(i);
			Variable_MCDsKrapidity = MCDsKrapidity->at(i);
			Variable_MCDsKmass = MCDsKmass->at(i);
			Variable_MCDsPieta = MCDsPieta->at(i);
			Variable_MCDsPiphi = MCDsPiphi->at(i);
			Variable_MCDsPipt = MCDsPipt->at(i);
			Variable_MCDsPienergy = MCDsPienergy->at(i);
			Variable_MCDsPip = MCDsPip->at(i);
			Variable_MCDsPiet = MCDsPiet->at(i);
			Variable_MCDsPirapidity = MCDsPirapidity->at(i);
			Variable_MCDsPimass = MCDsPimass->at(i);
	
			t_DsMC.Fill();	
		}

		//----------------------------------------
		//For For Matching
		//-----------------------------------------
		//double id = -1.;
		double deltaR =0.;
		std::vector<double> deltaR_vec;		
		//double minDeltaR = 0.;

		for(unsigned int k=0; k < MCDsmass->size(); k++)
		{
			for(unsigned int j=0; j < D0mass->size(); j++)
			{  //cout << "debug D* Matching 1 --------------------"<< endl;
													
				if( D0_VtxProb->at(j) < 0.01) continue;
				if( D0fromDSs3D->at(j) < SigCutD0fromDs) continue;
				if( Anglephi->at(j) < 0.99 ) continue;
				if ( D0pt->at(j) < 3. ) continue;
				if( (Dsmass->at(j) - D0mass->at(j)) > 0.16) continue; 

				//DeltaR Function
				cout << "MCDseta->at(i): " << MCDseta->at(k) << "----Dseta->at(j): " << Dseta->at(j);
				cout << "MCDsphi->at(i): " << MCDseta->at(k) << "----Dsphi->at(j): " << Dseta->at(j);

				deltaR = sqrt( pow( (MCDseta->at(k)-Dseta->at(j)) ,2) + pow( (MCDsphi->at(k)-Dsphi->at(j)) ,2));
				cout << "deltaR: " << deltaR;
				deltaR_vec.push_back(deltaR);
				//deltaR = sqrt( pow(MCDseta->at(j)-Dseta->at(j) ,2) + pow(MCDsphi->at(j)-Dsphi->at(j) ,2));
				cout << "debug D* Matching 2 --------------------"<< endl;			
			}	
		}
		if(deltaR_vec.size() > 0)cout << "deltaR_vec.size() --------------------"<< deltaR_vec.size() << endl;	
		//std::cout << "deltaR_vec[deltaR_vec] " << deltaR_vec[0] << endl;
		
		double minPos = 0, maxPos = 0;
    	for (unsigned i = 0; i < deltaR_vec.size(); ++i)
    	{ 	if (deltaR_vec[i] < deltaR_vec[minPos]) { minPos = i;} // Found a smaller min      
        	//if (deltaR_vec[i] > deltaR_vec[maxPos]) { maxPos = i;} // Found a bigger max
         //std::cout << " " << deltaR_vec[minPos] << " at position " << minPos << endl;    	
    	}
		//std::cout << "deltaR_vec[deltaR_vec] " << deltaR_vec[minPos] << endl;
    	
    	
		 if (debug)cout << "debug D* Matching process 6 --------------------"<< endl;
		//----------------------------------------
		//For D0 MC
		//-----------------------------------------	
		for(unsigned int h1=0; h1 < MCpromptD0eta->size(); h1++)
		{  
			//Variable_dScandsKpi = dScandsKpi->at(h1);
			Variable_MCpromptD0eta = MCpromptD0eta->at(h1);
			Variable_MCpromptD0phi = MCpromptD0phi->at(h1);
			Variable_MCpromptD0pt = MCpromptD0pt->at(h1);
			Variable_MCpromptD0energy = MCpromptD0energy->at(h1);
			Variable_MCpromptD0p = MCpromptD0p->at(h1);
			Variable_MCpromptD0et = MCpromptD0et->at(h1);
			Variable_MCpromptD0rapidity = MCpromptD0rapidity->at(h1);
			Variable_MCpromptD0mass = MCpromptD0mass->at(h1);
			Variable_MCpromptD0_Keta = MCpromptD0_Keta->at(h1);
			Variable_MCpromptD0_Kphi = MCpromptD0_Kphi->at(h1);
			Variable_MCpromptD0_Kpt = MCpromptD0_Kpt->at(h1);
			Variable_MCpromptD0_Kenergy = MCpromptD0_Kenergy->at(h1);
			Variable_MCpromptD0_Kp = MCpromptD0_Kp->at(h1);
			Variable_MCpromptD0_Ket = MCpromptD0_Ket->at(h1);
			Variable_MCpromptD0_Krapidity = MCpromptD0_Krapidity->at(h1);
			Variable_MCpromptD0_Kmass = MCpromptD0_Kmass->at(h1);
			Variable_MCpromptD0_Pieta = MCpromptD0_Pieta->at(h1);
			Variable_MCpromptD0_Piphi = MCpromptD0_Piphi->at(h1);
			Variable_MCpromptD0_Pipt = MCpromptD0_Pipt->at(h1);
			Variable_MCpromptD0_Pienergy = MCpromptD0_Pienergy->at(h1);
			Variable_MCpromptD0_Pip = MCpromptD0_Pip->at(h1);
			Variable_MCpromptD0_Piet = MCpromptD0_Piet->at(h1);
			Variable_MCpromptD0_Pirapidity = MCpromptD0_Pirapidity->at(h1);
			Variable_MCpromptD0_Pimass = MCpromptD0_Pimass->at(h1);
			Variable_MCpromptD0_DispAngle = MCpromptD0_DispAngle->at(h1);

			t_D0MC.Fill();	
		}

		//Mat D*
		/*for(unsigned int h2=0; h2 < DsetaMatching->size(); h2++)
		{  
			Variable_DsetaMatching = DsetaMatching->at(h2);
			Variable_MCDsetaMatching = MCDsetaMatching->at(h2);
			Variable_DsphiMatching = DsphiMatching->at(h2);
			Variable_MCDsphiMatching = MCDsphiMatching->at(h2);
			Variable_DsptMatching = DsptMatching->at(h2);
			Variable_MCDsptMatching = MCDsptMatching->at(h2); 
			Variable_D0fromDsmass = D0fromDsmass->at(h2);
			Variable_deltaRDsMatching = deltaRDsMatching->at(h2);

			if( DsptMatching->at(h2)>2 && MCDsptMatching->at(h2)>2 && abs(D0fromDsmass->at(h2) - 1.865)<0.024 )
			{				
				if( deltaRDsMatching->at(h2) < 0.01)
				{
					hRecEta_Rec1->Fill(abs(DsetaMatching->at(h2)));
					hRecEta_MC1->Fill(abs(MCDsetaMatching->at(h2)));
					hRecPhi_Rec1->Fill(DsphiMatching->at(h2));
					hRecPhi_MC1->Fill(MCDsphiMatching->at(h2));
					hRecPt_Rec1->Fill(DsptMatching->at(h2));
					hRecPt_MC1->Fill(MCDsptMatching->at(h2));																			
				}
				if( deltaRDsMatching->at(h2) < 0.02)
				{
					hRecEta_Rec2->Fill(abs(DsetaMatching->at(h2)));
					hRecEta_MC2->Fill(abs(MCDsetaMatching->at(h2)));
					hRecPhi_Rec2->Fill(DsphiMatching->at(h2));
					hRecPhi_MC2->Fill(MCDsphiMatching->at(h2));
					hRecPt_Rec2->Fill(DsptMatching->at(h2));
					hRecPt_MC2->Fill(MCDsptMatching->at(h2));	
				}
				if( deltaRDsMatching->at(h2) < 0.03 )
				{
					hRecEta_Rec3->Fill(abs(DsetaMatching->at(h2)));
					hRecEta_MC3->Fill(abs(MCDsetaMatching->at(h2)));
					hRecPhi_Rec3->Fill(DsphiMatching->at(h2));
					hRecPhi_MC3->Fill(MCDsphiMatching->at(h2));
					hRecPt_Rec3->Fill(DsptMatching->at(h2));
					hRecPt_MC3->Fill(MCDsptMatching->at(h2));
				}
				if( deltaRDsMatching->at(h2) < 0.04 )
				{
					hRecEta_Rec4->Fill(abs(DsetaMatching->at(h2)));
					hRecEta_MC4->Fill(abs(MCDsetaMatching->at(h2)));
					hRecPhi_Rec4->Fill(DsphiMatching->at(h2));
					hRecPhi_MC4->Fill(MCDsphiMatching->at(h2));
					hRecPt_Rec4->Fill(DsptMatching->at(h2));
					hRecPt_MC4->Fill(MCDsptMatching->at(h2));
				}
				if( deltaRDsMatching->at(h2) < 0.05 )
				{
					hRecEta_Rec5->Fill(abs(DsetaMatching->at(h2)));
					hRecEta_MC5->Fill(abs(MCDsetaMatching->at(h2)));
					hRecPhi_Rec5->Fill(DsphiMatching->at(h2));
					hRecPhi_MC5->Fill(MCDsphiMatching->at(h2));
					hRecPt_Rec5->Fill(DsptMatching->at(h2));
					hRecPt_MC5->Fill(MCDsptMatching->at(h2));
				}

				if( deltaRDsMatching->at(h2) < 0.06)
				{
					hRecEta_Rec6->Fill(abs(DsetaMatching->at(h2)));
					hRecEta_MC6->Fill(abs(MCDsetaMatching->at(h2)));
					hRecPhi_Rec6->Fill(DsphiMatching->at(h2));
					hRecPhi_MC6->Fill(MCDsphiMatching->at(h2));
					hRecPt_Rec6->Fill(DsptMatching->at(h2));
					hRecPt_MC6->Fill(MCDsptMatching->at(h2));
				}
			}
			t_DsMatching.Fill();	
		}

		for(unsigned int h3=0; h3 < D0etaMatching->size(); h3++)
		{  
			Variable_D0etaMatching = D0etaMatching->at(h3);
			Variable_MCD0etaMatching = MCD0etaMatching->at(h3);
			Variable_D0phiMatching = D0phiMatching->at(h3);
			Variable_MCD0phiMatching = MCD0phiMatching->at(h3);
			Variable_D0ptMatching = D0ptMatching->at(h3);
			Variable_MCD0ptMatching = MCD0ptMatching->at(h3);
			Variable_D0KtMatching = D0KtMatching->at(h3);
			Variable_D0SxyMatching = D0SxyMatching->at(h3);
			Variable_D0OpAngleMatching = D0OpAngleMatching->at(h3);
			Variable_deltaRD0Matching = deltaRD0Matching->at(h3);

			if(D0SxyMatching->at(h3)>2 && D0KtMatching->at(h3) > 0.3 && D0OpAngleMatching->at(h3) < 0.15)
			{
				if( deltaRD0Matching->at(h3) < 0.01)
				{
					hRecD0Eta_Rec1->Fill(abs(D0etaMatching->at(h3)));
					hRecD0Eta_MC1->Fill(abs(MCD0etaMatching->at(h3)));
					hRecD0Phi_Rec1->Fill(D0phiMatching->at(h3));
					hRecD0Phi_MC1->Fill(MCD0phiMatching->at(h3));
					hRecD0Pt_Rec1->Fill(D0ptMatching->at(h3));
					hRecD0Pt_MC1->Fill(MCD0ptMatching->at(h3));
				}
				if( deltaRD0Matching->at(h3) < 0.02)
				{
					hRecD0Eta_Rec2->Fill(abs(D0etaMatching->at(h3)));
					hRecD0Eta_MC2->Fill(abs(MCD0etaMatching->at(h3)));
					hRecD0Phi_Rec2->Fill(D0phiMatching->at(h3));
					hRecD0Phi_MC2->Fill(MCD0phiMatching->at(h3));
					hRecD0Pt_Rec2->Fill(D0ptMatching->at(h3));
					hRecD0Pt_MC2->Fill(MCD0ptMatching->at(h3));								
				}
				if( deltaRD0Matching->at(h3) < 0.03)
				{
					hRecD0Eta_Rec3->Fill(abs(D0etaMatching->at(h3)));
					hRecD0Eta_MC3->Fill(abs(MCD0etaMatching->at(h3)));
					hRecD0Phi_Rec3->Fill(D0phiMatching->at(h3));
					hRecD0Phi_MC3->Fill(MCD0phiMatching->at(h3));
					hRecD0Pt_Rec3->Fill(D0ptMatching->at(h3));
					hRecD0Pt_MC3->Fill(MCD0ptMatching->at(h3));                     
				}
				if( deltaRD0Matching->at(h3) < 0.04)
				{
					hRecD0Eta_Rec4->Fill(abs(D0etaMatching->at(h3)));
					hRecD0Eta_MC4->Fill(abs(MCD0etaMatching->at(h3)));
					hRecD0Phi_Rec4->Fill(D0phiMatching->at(h3));
					hRecD0Phi_MC4->Fill(MCD0phiMatching->at(h3));
					hRecD0Pt_Rec4->Fill(D0ptMatching->at(h3));
					hRecD0Pt_MC4->Fill(MCD0ptMatching->at(h3));              
				}
				if( deltaRD0Matching->at(h3) < 0.05)
				{
					hRecD0Eta_Rec5->Fill(abs(D0etaMatching->at(h3)));
					hRecD0Eta_MC5->Fill(abs(MCD0etaMatching->at(h3)));
					hRecD0Phi_Rec5->Fill(D0phiMatching->at(h3));
					hRecD0Phi_MC5->Fill(MCD0phiMatching->at(h3));
					hRecD0Pt_Rec5->Fill(D0ptMatching->at(h3));
					hRecD0Pt_MC5->Fill(MCD0ptMatching->at(h3));
				}

				if( deltaRD0Matching->at(h3) < 0.06)
				{
					hRecD0Eta_Rec6->Fill(abs(D0etaMatching->at(h3)));
					hRecD0Eta_MC6->Fill(abs(MCD0etaMatching->at(h3)));
					hRecD0Phi_Rec6->Fill(D0phiMatching->at(h3));
					hRecD0Phi_MC6->Fill(MCD0phiMatching->at(h3));
					hRecD0Pt_Rec6->Fill(D0ptMatching->at(h3));
					hRecD0Pt_MC6->Fill(MCD0ptMatching->at(h3));
				}

			}

			t_D0Matching.Fill();	
		}*/ 
			
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
	D0KpimassHisto->Write();
	t_analysis.Write();  //Write in the root file
	t_DstarWrongCombination.Write();
	t_D0analysis.Write();
	t_DsMC.Write();
	t_D0MC.Write();
	t_DsMatching.Write();
	t_D0Matching.Write();

	hRecEta_Rec1->Write();
	hRecEta_MC1->Write();
	hRecPhi_Rec1->Write();
	hRecPhi_MC1->Write();
	hRecPt_Rec1->Write();
	hRecPt_MC1->Write();

	hRecEta_Rec2->Write();
	hRecEta_MC2->Write();
	hRecPhi_Rec2->Write();
	hRecPhi_MC2->Write();
	hRecPt_Rec2->Write();
	hRecPt_MC2->Write();		

	hRecEta_Rec3->Write();
	hRecEta_MC3->Write();
	hRecPhi_Rec3->Write();
	hRecPhi_MC3->Write();
	hRecPt_Rec3->Write();
	hRecPt_MC3->Write();		

	hRecEta_Rec4->Write();
	hRecEta_MC4->Write();
	hRecPhi_Rec4->Write();
	hRecPhi_MC4->Write();
	hRecPt_Rec4->Write();
	hRecPt_MC4->Write();		

	hRecEta_Rec5->Write();
	hRecEta_MC5->Write();
	hRecPhi_Rec5->Write();
	hRecPhi_MC5->Write();
	hRecPt_Rec5->Write();
	hRecPt_MC5->Write();

	hRecEta_Rec6->Write();
	hRecEta_MC6->Write();
	hRecPhi_Rec6->Write();
	hRecPhi_MC6->Write();
	hRecPt_Rec6->Write();
	hRecPt_MC6->Write();

	//------------------------------------------
	hRecD0Eta_Rec1->Write();
	hRecD0Eta_MC1->Write();
	hRecD0Phi_Rec1->Write();
	hRecD0Phi_MC1->Write();
	hRecD0Pt_Rec1->Write();
	hRecD0Pt_MC1->Write();

	hRecD0Eta_Rec2->Write();
	hRecD0Eta_MC2->Write();
	hRecD0Phi_Rec2->Write();
	hRecD0Phi_MC2->Write();
	hRecD0Pt_Rec2->Write();
	hRecD0Pt_MC2->Write();

	hRecD0Eta_Rec3->Write();
	hRecD0Eta_MC3->Write();
	hRecD0Phi_Rec3->Write();
	hRecD0Phi_MC3->Write();
	hRecD0Pt_Rec3->Write();
	hRecD0Pt_MC3->Write();

	hRecD0Eta_Rec4->Write();
	hRecD0Eta_MC4->Write();
	hRecD0Phi_Rec4->Write();
	hRecD0Phi_MC4->Write();
	hRecD0Pt_Rec4->Write();
	hRecD0Pt_MC4->Write();

	hRecD0Eta_Rec5->Write();
	hRecD0Eta_MC5->Write();
	hRecD0Phi_Rec5->Write();
	hRecD0Phi_MC5->Write();
	hRecD0Pt_Rec5->Write();
	hRecD0Pt_MC5->Write();

	hRecD0Eta_Rec6->Write();
	hRecD0Eta_MC6->Write();
	hRecD0Phi_Rec6->Write();
	hRecD0Phi_MC6->Write();
	hRecD0Pt_Rec6->Write();
	hRecD0Pt_MC6->Write();

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout << "-------END PROGRAM-------------"<< endl;

}//end program
