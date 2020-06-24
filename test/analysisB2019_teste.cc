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
#include <TLatex.h>

//gSystem->AddIncludePath("/home/raphael/Desktop/analysis_2018");

using namespace std;
//--------------------------------------------------------------
//I F   I T   I S   M O N T E   C A R L O,  S E L E C T   Y E S
//--------------------------------------------------------------
#define GENvaraibles 1 // No=0 ; Yes=1 
//--------------------------------------------------------------
//S E C O N D A R Y   F U N C T I O N S
//--------------------------------------------------------------
TChain *callOneFile(const char* TITLES);
TChain *callTchain(const char* LOCAL);
TH1* makeTH1(const char* name, Double_t BIN, Double_t NMIN, Double_t NMAX, const char* TITLES, Color_t COLOR);
TH1* makeTH1Rec(const char* name, const int nbins_low , double* nbins_high );
//--------------------------------------------------------------
// M A I N   F U N C T I O N 
//--------------------------------------------------------------
void analysisB2019_OfflineCuts()
{

	clock_t tStart = clock();

	//Datasets: BparkingDataset1 - MC_DStarToD0Pi_D0KPi - MC_MinBias
	string Dataset = "MC_MinBias";

	string TriggerPath = ""; //Put "" for select no trigger //HLT_Mu9_IP6
	if (Dataset != "MC_DStarToD0Pi_D0KPi" and Dataset != "MC_MinBias"){TriggerPath = "HLT_Mu9_IP6";} //HLT_Mu9_IP6

	//Save the file
	string path = "/eos/user/r/ragomesd/crab/haddteste/";
	string path2 = "/eos/user/r/ragomesd/analysisB2019/canvas/";
	double SigCutD0fromDs, SigCutD0;

	if (Dataset != "MC_DStarToD0Pi_D0KPi"){SigCutD0fromDs = 3.; SigCutD0 = 5.; } //3
	if (Dataset == "MC_DStarToD0Pi_D0KPi" or Dataset == "MC_MinBias"){SigCutD0fromDs = 3.; SigCutD0 = 5.; }

	bool FlagDs = true;
	bool FlagDsWrong = true;
	bool FlagDsMC = true;
	bool FlagD0 = true;
	bool FlagD0MC = true;
	bool debug = false;
		
	//TChain* t1 = callTchain("/eos/user/r/ragomesd/crab/ParkingBPH2/ParkingBPH_Run2018A_MINIAOD/200105_042547/0000/");
	TChain* t1 = callTchain("./");

	//**********************************************************		
	Long64_t nentries = t1->GetEntries(); //Reading Number of tree entries of the file
	nentries = 100000; //Test
	cout<< "Number of tree entries: "<< nentries << endl;
	Long64_t partpercent = nentries*0.05; //Percent done to show
	double c = 299792458; //Light speed [m/s]
	
	//-------Reading the root file-------------------------------------	
	//TFile *f1 = new TFile("D0DstarDataBparking_7.root");
	//TTree *t1 = (TTree*)f1->Get("analysis/data");
	
	//-------Create root file-------------------------------------
	//TFile f_analysis("DMesonsHistograms.root","recreate"); //Creates root file
	TFile f_analysis(Form("%s%s.root", path.c_str(), Dataset.c_str() ),"recreate"); //Creates root file
	TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	TTree t_DstarWrongCombination("t_DstarWrongCombination","t_DstarWrongCombination");
	TTree t_D0analysis("t_D0analysis","t_D0analysis");
#if GENvaraibles == 0
#elif GENvaraibles == 1
	TTree t_DsMC("t_DsMC","t_DsMC"); 
	TTree t_D0MC("t_D0MC","t_D0MC");
	TTree t_DsMatching("t_DsMatching","t_DsMatching");
	TTree t_D0Matching("t_D0Matching","t_D0Matching");
#endif
	//Creates root file with Weighted Histograms
	TFile f_Puanalysis(Form("%sPU%s.root", path.c_str(), Dataset.c_str() ),"recreate"); 

	//-----------------------------------------------------
	//V A R I A B L E S
	//-----------------------------------------------------
	//It must be equal to "0", not "0.", even it is a double. Possible segmental error if you compile it with c++.
	vector<string>* NameOfFiredTriggers = 0;

	unsigned long 	NumberOfEvents, EventsAfterTrigger, PossibleDs, DsD0ProbVtx, DsD0Angle, D0fromDsMinusPDG, DsD0Sig, DsD0pt, DsD0deltaM, PossibleD0, D0ProbVtx, D0Angle, CountD0minusPDG, D0Sig, CountD0pt;
	NumberOfEvents = 0; EventsAfterTrigger = 0; PossibleDs = 0; DsD0ProbVtx = 0; DsD0Angle = 0; D0fromDsMinusPDG = 0; DsD0Sig = 0; DsD0pt = 0; DsD0deltaM = 0;	PossibleD0 = 0; D0ProbVtx = 0; D0Angle = 0; D0Sig = 0; CountD0minusPDG = 0; CountD0pt = 0;

	unsigned long GenDspT1, GenDspT2, GenDspT3, GenDspT4, GenDspT5, GenDspT6, GenDspT7, GenDspT8, GenDspT9, GenDsEta1, GenDsEta2, GenDsEta3, GenDsEta4, GenDsEta5, GenDsEta6, GenDsEta7, GenDsEta8, GenDsEta9, GenDsEta10;
	GenDspT1 = 0; GenDspT2 = 0; GenDspT3 = 0; GenDspT4 = 0; GenDspT5 = 0; GenDspT6 = 0; GenDspT7 = 0; GenDspT8 = 0; GenDspT9 = 0; GenDsEta1 = 0; GenDsEta2 = 0; GenDsEta3 = 0; GenDsEta4 = 0; GenDsEta5 = 0; GenDsEta6 = 0; GenDsEta7 = 0; GenDsEta8 = 0; GenDsEta9 = 0; GenDsEta10 = 0;

	unsigned long RecDspT1, RecDspT2, RecDspT3, RecDspT4, RecDspT5, RecDspT6, RecDspT7, RecDspT8, RecDspT9, RecDsEta1, RecDsEta2, RecDsEta3, RecDsEta4, RecDsEta5, RecDsEta6, RecDsEta7, RecDsEta8, RecDsEta9, RecDsEta10;
	RecDspT1 = 0; RecDspT2 = 0; RecDspT3 = 0; RecDspT4 = 0; RecDspT5 = 0; RecDspT6 = 0; RecDspT7 = 0; RecDspT8 = 0; RecDspT9 = 0; RecDsEta1 = 0; RecDsEta2 = 0; RecDsEta3 = 0; RecDsEta4 = 0; RecDsEta5 = 0; RecDsEta6 = 0; RecDsEta7 = 0; RecDsEta8 = 0; RecDsEta9 = 0; RecDsEta10 = 0;

	unsigned long GenD0pT1, GenD0pT2, GenD0pT3, GenD0pT4, GenD0pT5, GenD0pT6, GenD0pT7, GenD0pT8, GenD0pT9, GenD0Eta1, GenD0Eta2, GenD0Eta3, GenD0Eta4, GenD0Eta5, GenD0Eta6, GenD0Eta7, GenD0Eta8, GenD0Eta9, GenD0Eta10;
	GenD0pT1 = 0; GenD0pT2 = 0; GenD0pT3 = 0; GenD0pT4 = 0; GenD0pT5 = 0; GenD0pT6 = 0; GenD0pT7 = 0; GenD0pT8 = 0; GenD0pT9 = 0; GenD0Eta1 = 0; GenD0Eta2 = 0; GenD0Eta3 = 0; GenD0Eta4 = 0; GenD0Eta5 = 0; GenD0Eta6 = 0; GenD0Eta7 = 0; GenD0Eta8 = 0; GenD0Eta9 = 0; GenD0Eta10 = 0;

	unsigned long RecD0pT1, RecD0pT2, RecD0pT3, RecD0pT4, RecD0pT5, RecD0pT6, RecD0pT7, RecD0pT8, RecD0pT9, RecD0Eta1, RecD0Eta2, RecD0Eta3, RecD0Eta4, RecD0Eta5, RecD0Eta6, RecD0Eta7, RecD0Eta8, RecD0Eta9, RecD0Eta10;
	RecD0pT1 = 0; RecD0pT2 = 0; RecD0pT3 = 0; RecD0pT4 = 0; RecD0pT5 = 0; RecD0pT6 = 0; RecD0pT7 = 0; RecD0pT8 = 0; RecD0pT9 = 0; RecD0Eta1 = 0; RecD0Eta2 = 0; RecD0Eta3 = 0; RecD0Eta4 = 0; RecD0Eta5 = 0; RecD0Eta6 = 0; RecD0Eta7 = 0; RecD0Eta8 = 0;	RecD0Eta9 = 0; RecD0Eta10 = 0;
	//----------------------------------------
	//For D* 
	//-----------------------------------------
	vector<int>* FlagDstarfromB = 0;
	vector<double>* Dsmass = 0;
	vector<double>* Dseta = 0;
	vector<double>* Dsphi = 0;
	vector<double>* Dspt = 0;
	vector<double>* DSDeltaR = 0;

	vector<double>* D0mass = 0;
	vector<double>* D0eta = 0;
	vector<double>* D0phi = 0;
	vector<double>* D0pt = 0;
	vector<double>* D0fromDSdXY = 0; 
	vector<double>* D0fromDSsXY = 0;
	vector<double>* D0fromDSs3D = 0;
	vector<double>* D0fromDSd3D = 0;
	vector<double>* Anglephi = 0;
	vector<double>* D0_VtxProb = 0;

	/*vector<double>* TrkKpt = 0;
	vector<double>* TrkKnhits = 0;
	vector<double>* TrkKchi2 = 0;
	vector<double>* TrkKdxy = 0;
	vector<double>* TrkKdz = 0;
	vector<double>* TrkKeta = 0;
	vector<double>* TrkKphi = 0;	

	vector<double>* Trkpipt = 0;
	vector<double>* Trkpinhits = 0;
	vector<double>* Trkpichi2 = 0;
	vector<double>* Trkpidxy = 0;
	vector<double>* Trkpidz = 0;
	vector<double>* Trkpieta = 0;
	vector<double>* Trkpiphi = 0;

	vector<double>* TrkSpt = 0;
	vector<double>* TrkSnhits = 0;
	vector<double>* TrkSchi2 = 0;
	vector<double>* TrkSdxy = 0;
	vector<double>* TrkSdz = 0;
	vector<double>* TrkSeta = 0;
	vector<double>* TrkSphi = 0;*/
	//vector<double>* TrkScharge = 0;

	int DstarTotal = 0;
	int DstarSec = 0;
	double var_DsMinusD0 = 0.;
	double var_Dsmass = 0.;
	double var_Dseta = 0.;
	double var_Dsphi = 0.;
	double var_Dspt = 0.;

	double var_D0mass = 0.;
	double var_D0eta = 0.;
	double var_D0phi = 0.;
	double var_D0pt = 0.;
	double var_DSDeltaR = 0.;
	double var_D0fromDSsXY = 0.;
	double var_D0fromDSs3D = 0.;
	//double var_D0fromDSdXY = 0.; 
	//double var_D0fromDSd3D = 0.;
	double var_Anglephi = 0.;
	double var_D0_VtxProb = 0.;
	double var_Dslifetime = 0.;

	/*double var_TrkKpt = 0.;
	double var_TrkKnhits = 0.;
	double var_TrkKchi2 = 0.;
	double var_TrkKdxy = 0.;
	double var_TrkKdz = 0.;
	double var_TrkKeta = 0.;
	double var_TrkKphi = 0.;
	//double var_TrkScharge = 0.;

	double var_Trkpipt = 0.;
	double var_Trkpinhits = 0.;
	double var_Trkpichi2 = 0.;
	double var_Trkpidxy = 0.;
	double var_Trkpidz = 0.;
	double var_Trkpieta = 0.;
	double var_Trkpiphi = 0.;

	double var_TrkSpt = 0.;
	double var_TrkSnhits = 0.;
	double var_TrkSchi2 = 0.;
	double var_TrkSdxy = 0.;
	double var_TrkSdz = 0.;
	double var_TrkSeta = 0.;
	double var_TrkSphi = 0.;*/
	//----------------------------------------
	//For D* Wrong combination
	//-----------------------------------------
	vector<double>* D0massWrong = 0;
	vector<double>* D0etaWrong = 0;
	vector<double>* D0phiWrong = 0;
	vector<double>* D0ptWrong = 0;

	vector<double>* DsmassWrong = 0;
	vector<double>* DsetaWrong = 0;
	vector<double>* DsphiWrong = 0;
	vector<double>* DsptWrong = 0;

	//vector<double>* TrkScharge = 0;
	vector<double>* D0fromDSsXYWrong = 0;
	vector<double>* D0fromDSs3DWrong = 0;
	vector<double>* AnglephiWrong = 0;
	vector<double>* D0_VtxProbWrong = 0;

	double var_D0massWrong = 0.;
	double var_D0etaWrong = 0.;
	double var_D0phiWrong = 0.;
	double var_D0ptWrong = 0.;
	double var_DsMinusD0Wrong = 0.;
	double var_DsmassWrong = 0.;
	double var_DsetaWrong = 0.;
	double var_DsphiWrong = 0.;
	double var_DsptWrong = 0.;

	//double var_TrkSchargeWrong = 0.;
	double var_D0fromDSsXYWrong = 0.;
	double var_D0fromDSs3DWrong = 0.;
	double var_AnglephiWrong = 0.;
	double var_D0_VtxProbWrong = 0.;
	//----------------------------------------
	//For D0
	//-----------------------------------------
	vector<double>* D0Kpimass = 0;
	vector<double>* D0Kpi_VtxProb = 0;
	vector<double>* D0Kpis3D = 0;
	vector<double>* D0KpiDispAngle = 0;
	vector<double>* D0Kpipt = 0;
	vector<double>* D0Kpieta = 0;
	vector<double>* D0Kpiphi = 0;
	vector<double>* D0KpisXY = 0; 
	vector<double>* D0KpidXY = 0;
	vector<double>* D0Kpid3D = 0;
	vector<double>* D0Kpie3D = 0;
	vector<double>* D0_kT = 0;

	/*vector<double>* TrkD0Kdxy = 0;
	vector<double>* TrkD0Kdz = 0;
	vector<double>* TrkD0Kchi2 = 0;
	vector<double>* TrkD0Kpt = 0;
	vector<double>* TrkD0Keta = 0;
	vector<double>* TrkD0kphi = 0;
	vector<double>* TrkD0Knhits = 0;

	vector<double>* TrkD0pidxy = 0;
	vector<double>* TrkD0pidz = 0;
	vector<double>* TrkD0pichi2 = 0;
	vector<double>* TrkD0pipt = 0;
	vector<double>* TrkD0pieta = 0;
	vector<double>* TrkD0piphi = 0;
	vector<double>* TrkD0pinhits = 0;*/

	double var_D0Kpimass = 0.;
	double var_D0Kpi_VtxProb = 0.;
	double var_D0Kpis3D = 0.;
	double var_D0KpiDispAngle = 0.;

	double var_D0Kpipt = 0.;
	double var_D0Kpieta = 0.;
	double var_D0Kpiphi = 0.;
	double var_D0lifetime = 0.;
	double var_D0_kT = 0.;
	double var_D0KpidXY = 0.;
	double var_D0KpisXY = 0.;
	double var_D0Kpid3D = 0.;
	double var_D0Kpie3D = 0.;
	
	/*double var_TrkD0Kdxy = 0.;
	double var_TrkD0Kdz = 0.;
	double var_TrkD0Kchi2 = 0.;
	double var_TrkD0Kpt = 0.;
	double var_TrkD0Keta = 0.;
	double var_TrkD0kphi = 0.;
	double var_TrkD0Knhits = 0.;

	double var_TrkD0pidxy = 0.;
	double var_TrkD0pidz = 0.;
	double var_TrkD0pichi2 = 0.;
	double var_TrkD0pipt = 0.;
	double var_TrkD0pieta = 0.;
	double var_TrkD0piphi = 0.;
	double var_TrkD0pinhits = 0.;*/

#if GENvaraibles == 0	
#elif GENvaraibles == 1
	double PUWeight = 0.;
//======================================================
// MC Variables - D_star
//======================================================
	vector<double>* MCDseta = 0;
	vector<double>* MCDsphi = 0;
	vector<double>* MCDspt = 0;
	vector<double>* MCDsenergy = 0;
	vector<double>* MCDsp = 0;
	vector<double>* MCDset = 0;
	vector<double>* MCDsmass = 0;

	vector<double>* MCD0eta = 0;
	vector<double>* MCD0phi = 0;
	vector<double>* MCD0pt = 0;
	vector<double>* MCD0energy = 0;
	vector<double>* MCD0p = 0;
	vector<double>* MCD0et = 0;
	vector<double>* MCD0rapidity = 0;
	vector<double>* MCD0mass = 0;
	vector<double>* MCD0lifetime = 0;

	/*vector<double>* MCDsKphi = 0;
	vector<double>* MCDsKpt = 0;
	vector<double>* MCDsKenergy = 0;
	vector<double>* MCDsKp = 0;
	vector<double>* MCDsKet = 0;
	vector<double>* MCDsKrapidity = 0;
	vector<double>* MCDsKmass = 0;

	vector<double>* MCDsPieta = 0;
	vector<double>* MCDsPiphi = 0;
	vector<double>* MCDsPipt = 0;
	vector<double>* MCDsPienergy = 0;
	vector<double>* MCDsPip = 0;
	vector<double>* MCDsPiet = 0;
	vector<double>* MCDsPirapidity = 0;
	vector<double>* MCDsPimass = 0;*/

	/*double var_MCDseta = 0.;
	double var_MCDsphi = 0.;
	double var_MCDspt = 0.;
	double var_MCDsenergy = 0.;
	double var_MCDsp = 0.;
	double var_MCDset = 0.;
	double var_MCDsmass = 0.;

	double var_MCD0eta = 0.;
	double var_MCD0phi = 0.;
	double var_MCD0pt = 0.;
	double var_MCD0energy = 0.;
	double var_MCD0p = 0.;
	double var_MCD0et = 0.;
	double var_MCD0rapidity = 0.;
	double var_MCD0mass = 0.;*/
	double var_MCD0lifetime = 0.;

	/*double var_MCDsKphi = 0.;
	double var_MCDsKpt = 0.;
	double var_MCDsKenergy = 0.;
	double var_MCDsKp = 0.;
	double var_MCDsKet = 0.;
	double var_MCDsKrapidity = 0.;
	double var_MCDsKmass = 0.;

	double var_MCDsPieta = 0.;
	double var_MCDsPiphi = 0.;
	double var_MCDsPipt = 0.;
	double var_MCDsPienergy = 0.;
	double var_MCDsPip = 0.;
	double var_MCDsPiet = 0.;
	double var_MCDsPirapidity = 0.;
	double var_MCDsPimass = 0.;*/

//======================================================
// MC Variables - D0 
//======================================================
	vector<double>* MCpromptD0eta = 0;
	vector<double>* MCpromptD0phi = 0;
	vector<double>* MCpromptD0pt = 0;
	vector<double>* MCpromptD0energy = 0;
	vector<double>* MCpromptD0p = 0;
	vector<double>* MCpromptD0et = 0;
	vector<double>* MCpromptD0rapidity = 0;
	vector<double>* MCpromptD0mass = 0;
	vector<double>* MCpromptD0_DispAngle = 0;
	vector<double>* MCpromptD0lifetime = 0;
		
	/*vector<double>* MCpromptD0_Keta = 0;
	vector<double>* MCpromptD0_Kphi = 0;
	vector<double>* MCpromptD0_Kpt = 0;
	vector<double>* MCpromptD0_Kenergy = 0;
	vector<double>* MCpromptD0_Kp = 0;
	vector<double>* MCpromptD0_Ket = 0;
	vector<double>* MCpromptD0_Krapidity = 0;
	vector<double>* MCpromptD0_Kmass = 0;

	vector<double>* MCpromptD0_Pieta = 0;
	vector<double>* MCpromptD0_Piphi = 0;
	vector<double>* MCpromptD0_Pipt = 0;
	vector<double>* MCpromptD0_Pienergy = 0;
	vector<double>* MCpromptD0_Pip = 0;
	vector<double>* MCpromptD0_Piet = 0;
	vector<double>* MCpromptD0_Pirapidity = 0;
	vector<double>* MCpromptD0_Pimass = 0;*/

	/*double var_MCpromptD0eta = 0.;
	double var_MCpromptD0phi = 0.;
	double var_MCpromptD0pt = 0.;
	double var_MCpromptD0energy = 0.;
	double var_MCpromptD0p = 0.;
	double var_MCpromptD0et = 0.;
	double var_MCpromptD0rapidity = 0.;
	double var_MCpromptD0mass = 0.;
	double var_MCpromptD0_DispAngle = 0.;*/
	double var_MCpromptD0lifetime = 0.;

	/*double var_MCpromptD0_Keta = 0.;
	double var_MCpromptD0_Kphi = 0.;
	double var_MCpromptD0_Kpt = 0.;
	double var_MCpromptD0_Kenergy = 0.;
	double var_MCpromptD0_Kp = 0.;
	double var_MCpromptD0_Ket = 0.;
	double var_MCpromptD0_Krapidity = 0.;
	double var_MCpromptD0_Kmass = 0.;

	double var_MCpromptD0_Pieta = 0.;
	double var_MCpromptD0_Piphi = 0.;
	double var_MCpromptD0_Pipt = 0.;
	double var_MCpromptD0_Pienergy = 0.;
	double var_MCpromptD0_Pip = 0.;
	double var_MCpromptD0_Piet = 0.;
	double var_MCpromptD0_Pirapidity = 0.;
	double var_MCpromptD0_Pimass = 0.;*/
//======================================================
// MC Matching - D* 
//======================================================
	//vector<double>* MCD0lifetimeMatching = 0;
	/*vector<double>* DsetaMatching = 0;
	vector<double>* MCDsetaMatching = 0;
	vector<double>* DsphiMatching = 0;
	vector<double>* MCDsphiMatching = 0;
	vector<double>* DsptMatching = 0;
	vector<double>* MCDsptMatching = 0; 
	vector<double>* D0fromDsmass = 0;
	vector<double>* deltaRDsMatching = 0;*/

	double var_MCD0lifetimeMatching = 0.;
	/*double var_DsetaMatching = 0.;
	double var_MCDsetaMatching = 0.;
	double var_DsphiMatching = 0.;
	double var_MCDsphiMatching = 0.;
	double var_DsptMatching = 0.;
	double var_MCDsptMatching = 0.;
	double var_D0fromDsmass = 0.;*/
	double var_deltaRDsMatching = 0.;
	double DsdeltaEta = 0.;
	double DsdeltaPhi = 0.;
	double DsdeltaPt = 0.;
//======================================================
// MC Matching - D0 
//======================================================
	//vector<double>* MCpromptD0lifetimeMatching = 0;
	/*vector<double>* D0etaMatching = 0;
	vector<double>* MCD0etaMatching = 0;
	vector<double>* D0phiMatching = 0;
	vector<double>* MCD0phiMatching = 0;
	vector<double>* D0ptMatching = 0;
	vector<double>* MCD0ptMatching = 0;
	vector<double>* D0KtMatching = 0;
	vector<double>* D0SxyMatching = 0;
	vector<double>* D0OpAngleMatching = 0;
	vector<double>* deltaRD0Matching = 0;*/ 

	double var_MCpromptD0lifetimeMatching = 0.;
	/*double var_D0etaMatching = 0.;
	double var_MCD0etaMatching = 0.;
	double var_D0phiMatching = 0.;
	double var_MCD0phiMatching = 0.;
	double var_D0ptMatching = 0.;
	double var_MCD0ptMatching = 0.;
	double var_D0KtMatching = 0.;
	double var_D0SxyMatching = 0.;
	double var_D0OpAngleMatching = 0.;*/
	double var_deltaRD0Matching = 0.;
	double D0deltaEta = 0.;
	double D0deltaPhi = 0.;
	double D0deltaPt = 0.;
#endif

	string var_NameOfFiredTriggers; var_NameOfFiredTriggers.clear();
	
	//VARIABLES D*
	t_analysis.Branch("NameOfFiredTriggers",&var_NameOfFiredTriggers, "var_NameOfFiredTriggers/S");
	t_analysis.Branch("Dsmass",&var_Dsmass, "var_Dsmass/D");
	t_analysis.Branch("Dspt",&var_Dspt, "var_Dspt/D");
	t_analysis.Branch("Dseta",&var_Dseta, "var_Dseta/D");
	t_analysis.Branch("Dsphi",&var_Dsphi, "var_Dsphi/D");

	t_analysis.Branch("D0mass",&var_D0mass, "var_D0mass/D");
	t_analysis.Branch("DsMinusD0",&var_DsMinusD0, "var_DsMinusD0/D");
	t_analysis.Branch("D0pt",&var_D0pt, "var_D0pt/D");
	t_analysis.Branch("D0eta",&var_D0eta, "var_D0eta/D");
	t_analysis.Branch("D0phi",&var_D0phi, "var_D0phi/D");
	t_analysis.Branch("DSDeltaR",&var_DSDeltaR, "var_DSDeltaR/D");
	t_analysis.Branch("D0_VtxProb",&var_D0_VtxProb, "var_D0_VtxProb/D");
	t_analysis.Branch("Anglephi",&var_Anglephi, "var_Anglephi/D");
	t_analysis.Branch("D0fromDSsXY",&var_D0fromDSsXY, "var_D0fromDSsXY/D");
	t_analysis.Branch("D0fromDSs3D",&var_D0fromDSs3D, "var_D0fromDSs3D/D");
	t_analysis.Branch("Dslifetime",&var_Dslifetime, "var_Dslifetime/D");

	/*t_analysis.Branch("TrkKpt",&var_TrkKpt, "var_TrkKpt/D");
	t_analysis.Branch("TrkKnhits",&var_TrkKnhits, "var_TrkKnhits/D");
	t_analysis.Branch("TrkKchi2",&var_TrkKchi2, "var_TrkKchi2/D");
	t_analysis.Branch("TrkKdxy",&var_TrkKdxy, "var_TrkKdxy/D");
	t_analysis.Branch("TrkKdz",&var_TrkKdz, "var_TrkKdz/D");
	t_analysis.Branch("TrkKeta",&var_TrkKeta, "var_TrkKeta/D");      
	t_analysis.Branch("TrkKphi",&var_TrkKphi, "var_TrkKphi/D");

	t_analysis.Branch("Trkpipt",&var_Trkpipt, "var_Trkpipt/D");
	t_analysis.Branch("Trkpinhits",&var_Trkpinhits, "var_Trkpinhits/D");
	t_analysis.Branch("Trkpichi2",&var_Trkpichi2, "var_Trkpichi2/D");
	t_analysis.Branch("Trkpidxy",&var_Trkpidxy, "var_Trkpidxy/D");
	t_analysis.Branch("Trkpidz",&var_Trkpidz, "var_Trkpidz/D");
	t_analysis.Branch("Trkpieta",&var_Trkpieta, "var_Trkpieta/D");
	t_analysis.Branch("Trkpiphi",&var_Trkpiphi, "var_Trkpiphi/D");

	t_analysis.Branch("TrkSpt",&var_TrkSpt, "var_TrkSpt/D");
	t_analysis.Branch("TrkSnhits",&var_TrkSnhits, "var_TrkSnhits/D");
	t_analysis.Branch("TrkSchi2",&var_TrkSchi2, "var_TrkSchi2/D");
	t_analysis.Branch("TrkSdxy",&var_TrkSdxy, "var_TrkSdxy/D");	
	t_analysis.Branch("TrkSdz",&var_TrkSdz, "var_TrkSdz/D"); 
	t_analysis.Branch("TrkSeta",&var_TrkSeta, "var_TrkSeta/D");
	t_analysis.Branch("TrkSphi",&var_TrkSphi, "var_TrkSphi/D");*/
#if GENvaraibles == 0	
#elif GENvaraibles == 1
	t_analysis.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif
	//t_analysis.Branch("TrkScharge",&var_TrkScharge, "var_TrkScharge/D");

	//----------------------------------------
	//For D* Wrong combination
	//-----------------------------------------
	t_DstarWrongCombination.Branch("D0massWrong",&var_D0massWrong, "var_D0massWrong/D");
	t_DstarWrongCombination.Branch("D0ptWrong",&var_D0ptWrong, "var_D0ptWrong/D");
	t_DstarWrongCombination.Branch("D0etaWrong",&var_D0etaWrong, "var_D0etaWrong/D");
	t_DstarWrongCombination.Branch("D0phiWrong",&var_D0phiWrong, "var_D0phiWrong/D");

	t_DstarWrongCombination.Branch("DsMinusD0Wrong",&var_DsMinusD0Wrong, "var_DsMinusD0Wrong/D");
	
	t_DstarWrongCombination.Branch("DsmassWrong",&var_DsmassWrong, "var_DsmassWrong/D");
	t_DstarWrongCombination.Branch("DsptWrong",&var_DsptWrong, "var_DsptWrong/D");
	t_DstarWrongCombination.Branch("DsetaWrong",&var_DsetaWrong, "var_DsetaWrong/D");
	t_DstarWrongCombination.Branch("DsphiWrong",&var_DsphiWrong, "var_DsphiWrong/D");

	t_DstarWrongCombination.Branch("D0_VtxProbWrong",&var_D0_VtxProbWrong, "var_D0_VtxProbWrong/D");
	t_DstarWrongCombination.Branch("AnglephiWrong",&var_AnglephiWrong, "var_AnglephiWrong/D");
	t_DstarWrongCombination.Branch("D0fromDSsXYWrong",&var_D0fromDSsXYWrong, "var_D0fromDSsXYWrong/D");
	t_DstarWrongCombination.Branch("D0fromDSs3DWrong",&var_D0fromDSs3DWrong, "var_D0fromDSs3DWrong/D");
#if GENvaraibles == 0	
#elif GENvaraibles == 1
	t_DstarWrongCombination.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif

	//----------------------------------------
	//For D0 
	//-----------------------------------------
	t_D0analysis.Branch("D0Kpimass",&var_D0Kpimass, "var_D0Kpimass/D");
	t_D0analysis.Branch("D0Kpi_VtxProb",&var_D0Kpi_VtxProb, "var_D0Kpi_VtxProb/D");
	t_D0analysis.Branch("D0Kpis3D",&var_D0Kpis3D, "var_D0Kpis3D/D");
	t_D0analysis.Branch("D0KpiDispAngle",&var_D0KpiDispAngle, "var_D0KpiDispAngle/D");

	t_D0analysis.Branch("D0Kpipt",&var_D0Kpipt, "var_D0Kpipt/D");
	t_D0analysis.Branch("D0Kpieta",&var_D0Kpieta, "var_D0Kpieta/D");
	t_D0analysis.Branch("D0Kpiphi",&var_D0Kpiphi, "var_D0Kpiphi/D");
	t_D0analysis.Branch("D0lifetime",&var_D0lifetime, "var_D0lifetime/D");
	t_D0analysis.Branch("D0_kT",&var_D0_kT, "var_D0_kT/D");
	t_D0analysis.Branch("D0KpidXY",&var_D0KpidXY, "var_D0KpidXY/D"); 
	t_D0analysis.Branch("D0KpisXY",&var_D0KpisXY, "var_D0KpisXY/D");
	t_D0analysis.Branch("D0Kpid3D",&var_D0Kpid3D, "var_D0Kpid3D/D"); 
	t_D0analysis.Branch("D0Kpie3D",&var_D0Kpie3D, "var_D0Kpie3D/D");

	/*t_D0analysis.Branch("TrkD0Kdxy",&var_TrkD0Kdxy, "var_TrkD0Kdxy/D");
	t_D0analysis.Branch("TrkD0Kdz",&var_TrkD0Kdz, "var_TrkD0Kdz/D");
	t_D0analysis.Branch("TrkD0Kchi2",&var_TrkD0Kchi2, "var_TrkD0Kchi2/D");
	t_D0analysis.Branch("TrkD0Kpt",&var_TrkD0Kpt, "var_TrkD0Kpt/D");
	t_D0analysis.Branch("TrkD0Keta",&var_TrkD0Keta, "var_TrkD0Keta/D");
	t_D0analysis.Branch("TrkD0kphi",&var_TrkD0kphi, "var_TrkD0kphi/D");
	t_D0analysis.Branch("TrkD0Knhits",&var_TrkD0Knhits, "var_TrkD0Knhits/D");

	t_D0analysis.Branch("TrkD0pidxy",&var_TrkD0pidxy, "var_TrkD0pidxy/D");
	t_D0analysis.Branch("TrkD0pidz",&var_TrkD0pidz, "var_TrkD0pidz/D");
	t_D0analysis.Branch("TrkD0pichi2",&var_TrkD0pichi2, "var_TrkD0pichi2/D");
	t_D0analysis.Branch("TrkD0pipt",&var_TrkD0pipt, "var_TrkD0pipt/D");
	t_D0analysis.Branch("TrkD0pieta",&var_TrkD0pieta, "var_TrkD0pieta/D");
	t_D0analysis.Branch("TrkD0piphi",&var_TrkD0piphi, "var_TrkD0piphi/D");
	t_D0analysis.Branch("TrkD0pinhits",&var_TrkD0pinhits, "var_TrkD0pinhits/D");*/

#if GENvaraibles == 0	
#elif GENvaraibles == 1
	t_D0analysis.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif

#if GENvaraibles == 0
	
#elif GENvaraibles == 1
	//----------------------------------------
	//For D* MC
	//-----------------------------------------
	t_DsMC.Branch("MCD0lifetime",&var_MCD0lifetime, "var_MCD0lifetime/D");
	t_DsMC.Branch("PUWeight",&PUWeight, "PUWeight/D");

	/*t_DsMC.Branch("MCDseta",&var_MCDseta, "var_MCDseta/D");
	t_DsMC.Branch("MCDsphi",&var_MCDsphi, "var_MCDsphi/D");
	t_DsMC.Branch("MCDspt",&var_MCDspt, "var_MCDspt/D");
	t_DsMC.Branch("MCDsenergy",&var_MCDsenergy, "var_MCDsenergy/D");
	t_DsMC.Branch("MCDsp",&var_MCDsp, "var_MCDsp/D");
	t_DsMC.Branch("MCDset",&var_MCDset, "var_MCDset/D");
	t_DsMC.Branch("MCDsmass",&var_MCDsmass, "var_MCDsmass/D");

	t_DsMC.Branch("MCD0eta",&var_MCD0eta, "var_MCD0eta/D");
	t_DsMC.Branch("MCD0phi",&var_MCD0phi, "var_MCD0phi/D");
	t_DsMC.Branch("MCD0pt",&var_MCD0pt, "var_MCD0pt/D");
	t_DsMC.Branch("MCD0energy",&var_MCD0energy, "var_MCD0energy/D");
	t_DsMC.Branch("MCD0p",&var_MCD0p, "var_MCD0p/D");
	t_DsMC.Branch("MCD0et",&var_MCD0et, "var_MCD0et/D");
	t_DsMC.Branch("MCD0rapidity",&var_MCD0rapidity, "var_MCD0rapidity/D");
	t_DsMC.Branch("MCD0mass",&var_MCD0mass, "var_MCD0mass/D");*/

	/*t_DsMC.Branch("MCDsKphi",&var_MCDsKphi, "var_MCDsKphi/D");
	t_DsMC.Branch("MCDsKpt",&var_MCDsKpt, "var_MCDsKpt/D");
	t_DsMC.Branch("MCDsKenergy",&var_MCDsKenergy, "var_MCDsKenergy/D");
	t_DsMC.Branch("MCDsKp",&var_MCDsKp, "var_MCDsKp/D");
	t_DsMC.Branch("MCDsKet",&var_MCDsKet, "var_MCDsKet/D");
	t_DsMC.Branch("MCDsKrapidity",&var_MCDsKrapidity, "var_MCDsKrapidity/D");
	t_DsMC.Branch("MCDsKmass",&var_MCDsKmass, "var_MCDsKmass/D");

	t_DsMC.Branch("MCDsPieta",&var_MCDsPieta, "var_MCDsPieta/D");
	t_DsMC.Branch("MCDsPiphi",&var_MCDsPiphi, "var_MCDsPiphi/D");
	t_DsMC.Branch("MCDsPipt",&var_MCDsPipt, "var_MCDsPipt/D");
	t_DsMC.Branch("MCDsPienergy",&var_MCDsPienergy, "var_MCDsPienergy/D");
	t_DsMC.Branch("MCDsPip",&var_MCDsPip, "var_MCDsPip/D");
	t_DsMC.Branch("MCDsPiet",&var_MCDsPiet, "var_MCDsPiet/D");
	t_DsMC.Branch("MCDsPirapidity",&var_MCDsPirapidity, "var_MCDsPirapidity/D");
	t_DsMC.Branch("MCDsPimass",&var_MCDsPimass, "var_MCDsPimass/D");*/
	//----------------------------------------
	//For D0 MC
	//-----------------------------------------
	/*t_D0MC.Branch("MCpromptD0lifetime",&var_MCpromptD0lifetime, "var_MCpromptD0lifetime/D");
	t_D0MC.Branch("PUWeight",&PUWeight, "PUWeight/D");

	t_D0MC.Branch("MCpromptD0eta",&var_MCpromptD0eta, "var_MCpromptD0eta/D");
	t_D0MC.Branch("MCpromptD0phi",&var_MCpromptD0phi, "var_MCpromptD0phi/D");
	t_D0MC.Branch("MCpromptD0pt",&var_MCpromptD0pt, "var_MCpromptD0pt/D");
	t_D0MC.Branch("MCpromptD0energy",&var_MCpromptD0energy, "var_MCpromptD0energy/D");
	t_D0MC.Branch("MCpromptD0p",&var_MCpromptD0p, "var_MCpromptD0p/D");
	t_D0MC.Branch("MCpromptD0et",&var_MCpromptD0et, "var_MCpromptD0et/D");
	t_D0MC.Branch("MCpromptD0rapidity",&var_MCpromptD0rapidity, "var_MCpromptD0rapidity/D");
	t_D0MC.Branch("MCpromptD0mass",&var_MCpromptD0mass, "var_MCpromptD0mass/D");
	t_D0MC.Branch("MCpromptD0_DispAngle",&var_MCpromptD0_DispAngle, "var_MCpromptD0_DispAngle/D");*/


	/*t_D0MC.Branch("MCpromptD0_Keta",&var_MCpromptD0_Keta, "var_MCpromptD0_Keta/D");
	t_D0MC.Branch("MCpromptD0_Kphi",&var_MCpromptD0_Kphi, "var_MCpromptD0_Kphi/D");
	t_D0MC.Branch("MCpromptD0_Kpt",&var_MCpromptD0_Kpt, "var_MCpromptD0_Kpt/D");
	t_D0MC.Branch("MCpromptD0_Kenergy",&var_MCpromptD0_Kenergy, "var_MCpromptD0_Kenergy/D");
	t_D0MC.Branch("MCpromptD0_Kp",&var_MCpromptD0_Kp, "var_MCpromptD0_Kp/D");
	t_D0MC.Branch("MCpromptD0_Ket",&var_MCpromptD0_Ket, "var_MCpromptD0_Ket/D");
	t_D0MC.Branch("MCpromptD0_Krapidity",&var_MCpromptD0_Krapidity, "var_MCpromptD0_Krapidity/D");

	t_D0MC.Branch("MCpromptD0_Kmass",&var_MCpromptD0_Kmass, "var_MCpromptD0_Kmass/D");
	t_D0MC.Branch("MCpromptD0_Pieta",&var_MCpromptD0_Pieta, "var_MCpromptD0_Pieta/D");
	t_D0MC.Branch("MCpromptD0_Piphi",&var_MCpromptD0_Piphi, "var_MCpromptD0_Piphi/D");
	t_D0MC.Branch("MCpromptD0_Pipt",&var_MCpromptD0_Pipt, "var_MCpromptD0_Pipt/D");
	t_D0MC.Branch("MCpromptD0_Pienergy",&var_MCpromptD0_Pienergy, "var_MCpromptD0_Pienergy/D");
	t_D0MC.Branch("MCpromptD0_Pip",&var_MCpromptD0_Pip, "var_MCpromptD0_Pip/D");
	t_D0MC.Branch("MCpromptD0_Piet",&var_MCpromptD0_Piet, "var_MCpromptD0_Piet/D");
	t_D0MC.Branch("MCpromptD0_Pirapidity",&var_MCpromptD0_Pirapidity, "var_MCpromptD0_Pirapidity/D");
	t_D0MC.Branch("MCpromptD0_Pimass",&var_MCpromptD0_Pimass, "var_MCpromptD0_Pimass/D");*/
	//----------------------------------------
	//For D* Matching
	//-----------------------------------------
	t_DsMatching.Branch("MCD0lifetimeMatching",&var_MCD0lifetimeMatching, "var_MCD0lifetimeMatching/D");
	t_DsMatching.Branch("PUWeight",&PUWeight, "PUWeight/D");
	/*t_DsMatching.Branch("DsetaMatching",&var_DsetaMatching, "var_DsetaMatching/D");
	t_DsMatching.Branch("MCDsetaMatching",&var_MCDsetaMatching, "var_MCDsetaMatching/D");
	t_DsMatching.Branch("DsphiMatching",&var_DsphiMatching, "var_DsphiMatching/D");
	t_DsMatching.Branch("MCDsphiMatching",&var_MCDsphiMatching, "var_MCDsphiMatching/D");
	t_DsMatching.Branch("DsptMatching",&var_DsptMatching, "var_DsptMatching/D");
	t_DsMatching.Branch("MCDsptMatching",&var_MCDsptMatching, "var_MCDsptMatching/D");
	t_DsMatching.Branch("D0fromDsmass",&var_D0fromDsmass, "var_D0fromDsmass/D");
	t_DsMatching.Branch("deltaRDsMatching",&var_deltaRDsMatching, "var_deltaRDsMatching/D");
	t_DsMatching.Branch("DsdeltaEta",&DsdeltaEta, "DsdeltaEta/D");
	t_DsMatching.Branch("DsdeltaPhi",&DsdeltaPhi, "DsdeltaPhi/D"); 
	t_DsMatching.Branch("DsdeltaPt",&DsdeltaPt, "DsdeltaPt/D");*/
	//----------------------------------------
	//For D0 Matching
	//----------------------------------------- 
	t_D0Matching.Branch("MCpromptD0lifetimeMatching",&var_MCpromptD0lifetimeMatching, "var_MCpromptD0lifetimeMatching/D");
	t_D0Matching.Branch("PUWeight",&PUWeight, "PUWeight/D");
	/*t_D0Matching.Branch("D0etaMatching",&var_D0etaMatching, "var_D0etaMatching/D");
	t_D0Matching.Branch("MCD0etaMatching",&var_MCD0etaMatching, "var_MCD0etaMatching/D");
	t_D0Matching.Branch("D0phiMatching",&var_D0phiMatching, "var_D0phiMatching/D");
	t_D0Matching.Branch("MCD0phiMatching",&var_MCD0phiMatching, "var_MCD0phiMatching/D");
	t_D0Matching.Branch("D0ptMatching",&var_D0ptMatching, "var_D0ptMatching/D");
	t_D0Matching.Branch("MCD0ptMatching",&var_MCD0ptMatching, "var_MCD0ptMatching/D");
	t_D0Matching.Branch("D0KtMatching",&var_D0KtMatching, "var_D0KtMatching/D");
	t_D0Matching.Branch("D0SxyMatching",&var_D0SxyMatching, "var_D0SxyMatching/D");
	t_D0Matching.Branch("D0OpAngleMatching",&var_D0OpAngleMatching, "var_D0OpAngleMatching/D");
	t_D0Matching.Branch("deltaRD0Matching",&var_deltaRD0Matching, "var_deltaRD0Matching/D");
	t_D0Matching.Branch("D0deltaEta",&D0deltaEta, "D0deltaEta/D");
	t_D0Matching.Branch("D0deltaPhi",&D0deltaPhi, "D0deltaPhi/D"); 
	t_D0Matching.Branch("D0deltaPt",&D0deltaPt, "D0deltaPt/D");*/
#endif
	if (debug)cout << "debug 2 --------------------" << endl;
	//--------------------------------------------------
	//Creating Histgrams
	//---------------------------------------------------
	//TH1* makeTH1(const char* name, Double_t BIN, Double_t NMIN, Double_t NMAX, const char* TITLES, Color_t COLOR,const char* TYPE)
	TH1 *D0massHisto = makeTH1("D0massHisto", 100, 1.76, 1.96, "Invariant Mass of the D0(From D*) ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassHisto = makeTH1("DsmassHisto", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *D0fromDSs3DHisto = makeTH1("D0fromDSs3DHisto", 100,0,5, "Significance of the D0(From D*) ; Significance ; Events ", kRed);
	TH1 *TotalD0fromDSs3DHisto = makeTH1("TotalD0fromDSs3DHisto", 100,0,5, "Significance of the D0(From D*) ; Significance ; Events ", kRed);
	TH1 *DsMinusD0Histo = makeTH1("DsMinusD0Histo", 100,0.14,0.16, "#Delta m = m(K#pi#pi_{slow}) - m(K#pi) ; #Delta m ; Events", kRed);
	TH1 *D0VtxProbHisto = makeTH1("D0VtxProbHisto", 100,0,1, "D0(From D*) Vtx Prob ; Probability ; Events ", kRed);
	TH1 *DslifetimeHisto = makeTH1("DslifetimeHisto", 100,0,2*(pow(10,-12)), "D0(From D*) lifetime ; time [s] ; Events ", kRed);

	//Study of Mass with different D0 significance cuts
	TH1 *DsmassHisto0 = makeTH1("DsmassHisto0", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassHisto1 = makeTH1("DsmassHisto1", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassHisto2 = makeTH1("DsmassHisto2", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassHisto3 = makeTH1("DsmassHisto3", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassHisto4 = makeTH1("DsmassHisto4", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassHisto5 = makeTH1("DsmassHisto5", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassHisto6 = makeTH1("DsmassHisto6", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);

	TH1 *DsmassPtHisto1 = makeTH1("DsmassPtHisto1", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassPtHisto2 = makeTH1("DsmassPtHisto2", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassPtHisto3 = makeTH1("DsmassPtHisto3", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassPtHisto4 = makeTH1("DsmassPtHisto4", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassPtHisto5 = makeTH1("DsmassPtHisto5", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassPtHisto6 = makeTH1("DsmassPtHisto6", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassPtHisto7 = makeTH1("DsmassPtHisto7", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassPtHisto8 = makeTH1("DsmassPtHisto8", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);
	TH1 *DsmassPtHisto9 = makeTH1("DsmassPtHisto9", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV] ; Events ", kRed);

	TH1 *D0KpimassHisto = makeTH1("D0KpimassHisto", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0Kpi_VtxProbHisto = makeTH1("D0Kpi_VtxProbHisto", 100,0,1, "D0Kpi_VtxProbHisto D0 ; Probability ; Events ", kRed);
	TH1 *TotalD0Kpis3DHistoHisto = makeTH1("TotalD0Kpis3DHistoHisto", 100,0,5, "Significance of the D0 ; Significance ; Events ", kRed);
	TH1 *D0Kpis3DHisto = makeTH1("D0Kpis3DHisto", 100,0,6, "Significance of the D0 ; Significance ; Events ", kRed);
	TH1 *D0lifetimeHisto = makeTH1("D0lifetimeHisto", 100,0,2*(pow(10,-12)), " D0 lifetime ; time [s] ; Events ", kRed);

	//Study of Mass with different D0 significance cuts
	TH1 *D0KpimassHisto0 = makeTH1("D0KpimassHisto0", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto1 = makeTH1("D0KpimassHisto1", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto2 = makeTH1("D0KpimassHisto2", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto3 = makeTH1("D0KpimassHisto3", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto4 = makeTH1("D0KpimassHisto4", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto5 = makeTH1("D0KpimassHisto5", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto6 = makeTH1("D0KpimassHisto6", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto7 = makeTH1("D0KpimassHisto7", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto8 = makeTH1("D0KpimassHisto8", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto9 = makeTH1("D0KpimassHisto9", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassHisto10 = makeTH1("D0KpimassHisto10", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);

	TH1 *D0KpimassPtHisto1 = makeTH1("D0KpimassPtHisto1", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassPtHisto2 = makeTH1("D0KpimassPtHisto2", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassPtHisto3 = makeTH1("D0KpimassPtHisto3", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassPtHisto4 = makeTH1("D0KpimassPtHisto4", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassPtHisto5 = makeTH1("D0KpimassPtHisto5", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassPtHisto6 = makeTH1("D0KpimassPtHisto6", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassPtHisto7 = makeTH1("D0KpimassPtHisto7", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassPtHisto8 = makeTH1("D0KpimassPtHisto8", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *D0KpimassPtHisto9 = makeTH1("D0KpimassPtHisto9", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);

	//D* PU weight
	TH1 *hDsmass = makeTH1("hDsmass", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV/c2] ; Events ", kRed);
	TH1 *hDseta = makeTH1("hDseta", 100,-4,4, "Pseudo-rapidity distribuition of the D0(From D*) ; #eta ; Events", kRed);
	TH1 *hDsphi = makeTH1("hDsphi", 100,-4,4, "#Phi distribuition of the D* ; #Phi ; Events ", kRed);
	TH1 *hDspt = makeTH1("hDspt", 100,0,30, "pT distribuition of the D* ; pT [GeV] ; Events", kRed);
	TH1 *hDslifetime = makeTH1("hDslifetime",100,0,2*(pow(10,-12)), "lifetime of the D0(FromD*) ; lifetime [s] ; Events ", kRed);

	TH1 *hDsMinusD0 = makeTH1("hDsMinusD0", 100,0.13,0.17, "#Delta m = m(K#pi#pi_{slow}) - m(K#pi) ; #Delta m [GeV/c2]; Events", kRed);
	TH1 *hD0mass = makeTH1("hD0mass", 100,1.76,1.96, "Invariant Mass of the D0(FromD*) ; Mass [GeV] ; Events ", kRed);
	TH1 *hD0eta = makeTH1("hD0eta", 100,-4,4, "Invariant Mass of the D0(FromD*) ; Mass [GeV] ; Events ", kRed);
	TH1 *hD0phi = makeTH1("hD0phi", 100,-4,4, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *hD0pt = makeTH1("hD0pt", 100,0,30, "Invariant Mass of the D0(FromD*) ; Mass [GeV] ; Events ", kRed);

	TH1 *hTrkKpt = makeTH1("hTrkKpt", 100,0,10, "pt distribuition of the Kaon ; pT [GeV] ; Events ", kRed);
	TH1 *hTrkKchi2 = makeTH1("hTrkKchi2", 100,0,4, "Chi2 distribuition of the Kaon ; #Chi^{2} ; Events", kRed);
	TH1 *hTrkKnhits = makeTH1("hTrkKnhits", 100,0,50, "nhits distribuition of the Pion ; Number of hits ; Events", kRed);
	TH1 *hTrkKdz = makeTH1("hTrkKdz", 100,-0.1,0.1, "dz distribuition of the Kaon ; dz [cm] ; Events", kRed);
	TH1 *hTrkKdxy = makeTH1("hTrkKdxy", 100,-0.1,0.1, "dxy distribuition of the Kaon ; dz [cm] ; Events", kRed);
	TH1 *hTrkKeta = makeTH1("hTrkKeta", 100,-4,4, "eta distribuition of the Kaon ; #eta ; Events", kRed);
	TH1 *hTrkKphi = makeTH1("hTrkKphi", 100,-4,4, "phi distribuition of the Kaon ; #Phi ; Events", kRed);

	TH1 *hTrkpipt = makeTH1("hTrkpipt", 100,0,10, "pt distribuition of the Pion ; pT [GeV] ; Events ", kRed);
	TH1 *hTrkpichi2 = makeTH1("hTrkpichi2", 100,0,4, "Chi2 distribuition of the Pion ; #Chi^{2} ; Events", kRed);
	TH1 *hTrkpinhits = makeTH1("hTrkpinhits", 100,0,50, "nhits distribuition of the Pion ; Number of hits ; Events", kRed);
	TH1 *hTrkpidz = makeTH1("hTrkpidz", 100,-0.1,0.1, "dz distribuition of the Pion ; dz [cm] ; Events", kRed);
	TH1 *hTrkpidxy = makeTH1("hTrkpidxy", 100,-0.1,0.1, "dxy distribuition of the Pion ; dz [cm] ; Events", kRed);
	TH1 *hTrkpieta = makeTH1("hTrkpieta", 100,-4,4, "eta distribuition of the Pion ; #eta ; Events", kRed);
	TH1 *hTrkpiphi = makeTH1("hTrkpiphi", 100,-4,4, "phi distribuition of the Pion ; #Phi ; Events", kRed);

	TH1 *hTrkSpt = makeTH1("hTrkSpt", 100,0,10, "pt distribuition of the SlowPion ; pT [GeV] ; Events ", kRed);
	TH1 *hTrkSchi2 = makeTH1("hTrkSchi2", 100,0,4, "Chi2 distribuition of the SlowPion ; #Chi^{2} ; Events", kRed);
	TH1 *hTrkSnhits = makeTH1("hTrkSnhits", 100,0,50, "nhits distribuition of the SlowPion ; Number of hits ; Events", kRed);
	TH1 *hTrkSdz = makeTH1("hTrkSdz", 100,-0.1,0.1, "dz distribuition of the SlowPion ; dz [cm] ; Events", kRed);
	TH1 *hTrkSdxy = makeTH1("hTrkSdxy", 100,-0.1,0.1, "dxy distribuition of the SlowPion ; dz [cm] ; Events", kRed);
	TH1 *hTrkSeta = makeTH1("hTrkSeta", 100,-4,4, "eta distribuition of the SlowPion ; #eta ; Events", kRed);
	TH1 *hTrkSphi = makeTH1("hTrkSphi", 100,-4,4, "phi distribuition of the SlowPion ; #Phi ; Events", kRed);

	//D0 PU weight
	TH1 *hD0Kpimass = makeTH1("hD0Kpimass", 100,1.76,1.96, "Invariant Mass of the From D0 ; Mass [GeV/c2] ; Events ", kRed);
	TH1 *hD0Kpieta = makeTH1("hD0Kpieta", 100,-4,4, "Pseudo-rapidity distribuition of the D0(From D*) ; #eta ; Events", kRed);
	TH1 *hD0Kpiphi = makeTH1("hD0Kpiphi", 100,-4,4, "#Phi distribuition of the D0 ; #Phi ; Events ", kRed);
	TH1 *hD0Kpipt = makeTH1("hD0Kpipt", 100,0,30, "pT distribuition of the D0 ; pT [GeV] ; Events", kRed);
	TH1 *hD0lifetime = makeTH1("hD0lifetime",  100,0,2*(pow(10,-12)), "lifetime of the D0 ; lifetime [s] ; Events ", kRed);
	
	TH1 *hTrkD0Kpt = makeTH1("hTrkD0Kpt", 100,0,10, "pt distribuition of the Kaon ; pT [GeV] ; Events ", kRed);
	TH1 *hTrkD0Kchi2 = makeTH1("hTrkD0Kchi2", 100,0,4, "Chi2 distribuition of the Kaon ; #Chi^{2} ; Events", kRed);
	TH1 *hTrkD0Knhits = makeTH1("hTrkD0Knhits", 100,0,50, "nhits distribuition of the Pion ; Number of hits ; Events", kRed);
	TH1 *hTrkD0Kdz = makeTH1("hTrkD0Kdz", 100,-0.1,0.1, "dz distribuition of the Kaon ; dz [cm] ; Events", kRed);
	TH1 *hTrkD0Kdxy = makeTH1("hTrkD0Kdxy", 100,-0.1,0.1, "dxy distribuition of the Kaon ; dz [cm] ; Events", kRed);
	TH1 *hTrkD0Keta = makeTH1("hTrkD0Keta", 100,-4,4, "eta distribuition of the Kaon ; #eta ; Events", kRed);
	TH1 *hTrkD0kphi = makeTH1("hTrkD0kphi", 100,-4,4, "phi distribuition of the Kaon ; #Phi ; Events", kRed);

	TH1 *hTrkD0pipt = makeTH1("hTrkD0pipt", 100,0,10, "pt distribuition of the Pion ; pT [GeV] ; Events ", kRed);
	TH1 *hTrkD0pichi2 = makeTH1("hTrkD0pichi2", 100,0,4, "Chi2 distribuition of the Pion ; #Chi^{2} ; Events", kRed);
	TH1 *hTrkD0pinhits = makeTH1("hTrkD0pinhits", 100,0,50, "nhits distribuition of the Pion ; Number of hits ; Events", kRed);
	TH1 *hTrkD0pidz = makeTH1("hTrkD0pidz", 100,-0.1,0.1, "dz distribuition of the Pion ; dz [cm] ; Events", kRed);
	TH1 *hTrkD0pidxy = makeTH1("hTrkD0pidxy", 100,-0.1,0.1, "dxy distribuition of the Pion ; dz [cm] ; Events", kRed);
	TH1 *hTrkD0pieta = makeTH1("hTrkD0pieta", 100,-4,4, "eta distribuition of the Pion ; #eta ; Events", kRed);
	TH1 *hTrkD0piphi = makeTH1("hTrkD0piphi", 100,-4,4, "phi distribuition of the Pion ; #Phi ; Events", kRed);
	
	//TH1 *D0KpiDispAngleHisto = makeTH1("D0KpiDispAngleHisto", 100,0,1, "Invariant Mass of the D0(Prompt) ; Mass [GeV] ; Events ", kRed);

	/*TH1 *hRecEta_MC = makeTH1("hRecEta_MC", 10,0,2.5, "hRecEta_MC ; #eta ; Events ", kRed);
	TH1 *hRecPhi_MC = makeTH1("hRecPhi_MC", 10,0,5, "hRecPhi_MC ; #phi ; Events ", kRed);
	//TH1 *hRecPt_MC = makeTH1("hRecPt_MC", 10,0,40, "hRecPt_MC ; pT [GeV] ; Events ", kRed);
	
	TH1 *hRecEta_Rec1 = makeTH1("hRecEta_Rec1", 10,0,2.5, "hRecEta_Rec1 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec1 = makeTH1("hRecPhi_Rec1", 10,0,5, "hRecPhi_Rec1 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec1 = makeTH1("hRecPt_Rec1", 10,0,40, "hRecPt_Rec1 ; pT [GeV] ; Events ", kRed);

	TH1 *hRecEta_Rec2 = makeTH1("hRecEta_Rec2", 10,0,2.5, "hRecEta_Rec2 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec2 = makeTH1("hRecPhi_Rec2", 10,0,5, "hRecPhi_Rec2 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec2 = makeTH1("hRecPt_Rec2", 10,0,40, "hRecPt_Rec2 ; pT [GeV] ; Events ", kRed);

	TH1 *hRecEta_Rec3 = makeTH1("hRecEta_Rec3", 10,0,2.5, "hRecEta_Rec3 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec3 = makeTH1("hRecPhi_Rec3", 10,0,5, "hRecPhi_Rec3 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec3 = makeTH1("hRecPt_Rec3", 10,0,40, "hRecPt_Rec3 ; pT [GeV] ; Events ", kRed);

	TH1 *hRecEta_Rec4 = makeTH1("hRecEta_Rec4", 10,0,2.5, "hRecEta_Rec4 ; #eta ; Events ", kRed);
	TH1 *hRecPhi_Rec4 = makeTH1("hRecPhi_Rec4", 10,0,5, "hRecPhi_Rec4 ; #phi ; Events ", kRed);
	TH1 *hRecPt_Rec4 = makeTH1("hRecPt_Rec4", 10,0,40, "hRecPt_Rec4 ; pT [GeV]; Events ", kRed);*/

#if GENvaraibles == 0	
#elif GENvaraibles == 1
	TH2D* hdeltaRDstar = new TH2D("hdeltaRDstar", "hdeltaRDstar", 100, 1.94, 2.1, 100, 1.94, 2.1) ;
	hdeltaRDstar->SetTitle("hdeltaRDstar ; Rec Mass[GeV]; Gen Mass[GeV]"); hdeltaRDstar->SetMarkerStyle(1);
	hdeltaRDstar->Sumw2();	//hdeltaRDstar->SetFillColor(COLOR);

	TH2D* hdeltaRDstarPt = new TH2D("hdeltaRDstarPt", "hdeltaRDstarPt", 100, -1, 40., 100, -1, 40. ) ;
	hdeltaRDstarPt->SetTitle("hdeltaRDstarPt ; Rec pT[GeV]; Gen pT[GeV]"); hdeltaRDstarPt->SetMarkerStyle(1);
	hdeltaRDstarPt->Sumw2();	//hdeltaRDstarPt->SetFillColor(COLOR);

	TH2D* hdeltaRDstarEta = new TH2D("hdeltaRDstarEta", "hdeltaRDstarEta", 100, -1, 3., 100, -1, 3.) ;
	hdeltaRDstarEta->SetTitle("hdeltaRDstarEta ; Rec #eta; Gen #eta"); hdeltaRDstarEta->SetMarkerStyle(1);
	hdeltaRDstarEta->Sumw2();	//hdeltaRDstarEta->SetFillColor(COLOR);

	TH2D* hdeltaRDstarTime = new TH2D("hdeltaRDstarTime", "hdeltaRDstarTime", 100, 0, 1*pow(10,-12), 100, 0, 1*pow(10,-12)) ;
	hdeltaRDstarTime->SetTitle("hdeltaRDstarTime ; Rec t[s]; Gen t[s]"); hdeltaRDstarTime->SetMarkerStyle(1);
	hdeltaRDstarTime->Sumw2();	//hdeltaRDstarTime->SetFillColor(COLOR);

	TH2D* hdeltaRD0 = new TH2D("hdeltaRD0", "hdeltaRD0", 100, 1.76, 1.96, 100, 1.76, 1.96) ;
	hdeltaRD0->SetTitle("hdeltaRD0 ; Rec Mass[GeV]; Gen Mass[GeV]"); hdeltaRD0->SetMarkerStyle(1);
	hdeltaRD0->Sumw2();	//hdeltaRD0->SetFillColor(COLOR);

	TH2D* hdeltaRD0Pt = new TH2D("hdeltaRD0Pt", "hdeltaRD0Pt", 100, -1, 30., 100, -1, 30.) ;
	hdeltaRD0Pt->SetTitle("hdeltaRD0Pt ; Rec pT[GeV]; Gen pT[GeV]"); hdeltaRD0Pt->SetMarkerStyle(1);
	hdeltaRD0Pt->Sumw2();	//hdeltaRD0Pt->SetFillColor(COLOR); 

	TH2D* hdeltaRD0Eta = new TH2D("hdeltaRD0Eta", "hdeltaRD0Eta", 100, -1, 3., 100, -1, 3.) ;
	hdeltaRD0Eta->SetTitle("hdeltaRD0Eta ; Rec #eta; Gen #eta"); hdeltaRD0Eta->SetMarkerStyle(1);
	hdeltaRD0Eta->Sumw2();	//hdeltaRD0Eta->SetFillColor(COLOR);

	TH2D* hdeltaRD0Time = new TH2D("hdeltaRD0Time", "hdeltaRD0Time", 100, 0, 3*pow(10,-12), 100, 0, 3*pow(10,-12)) ;
	hdeltaRD0Time->SetTitle("hdeltaRD0Time ; Rec t[s]; Gen t[s]"); hdeltaRD0Time->SetMarkerStyle(1);
	hdeltaRD0Time->Sumw2();	//hdeltaRD0Time->SetFillColor(COLOR);

	const int nbinsPt_left = 9;
	double binsPt_left[nbinsPt_left+1] = {4., 5., 6., 7., 8., 12., 16., 24., 40., 100.};
	const int nbinsEta_left = 10;
	double binsEta_left[nbinsEta_left+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1};

	//For D*
	TH1 *hRecEta_MC = makeTH1Rec("hRecEta_MC", nbinsEta_left, binsEta_left );
	TH1 *hRecPt_MC = makeTH1Rec("hRecPt_MC", nbinsPt_left, binsPt_left );

	TH1 *hRecEta_Rec1 = makeTH1Rec("hRecEta_Rec1", nbinsEta_left, binsEta_left );
	TH1 *hRecPt_Rec1 = makeTH1Rec("hRecPt_Rec1", nbinsPt_left, binsPt_left );

	TH1 *hRecEta_Rec2 = makeTH1Rec("hRecEta_Rec2", nbinsEta_left, binsEta_left );
	TH1 *hRecPt_Rec2 = makeTH1Rec("hRecPt_Rec2", nbinsPt_left, binsPt_left );

	TH1 *hRecEta_Rec3 = makeTH1Rec("hRecEta_Rec3", nbinsEta_left, binsEta_left );
	TH1 *hRecPt_Rec3 = makeTH1Rec("hRecPt_Rec3", nbinsPt_left, binsPt_left );

	TH1 *hRecEta_Rec4 = makeTH1Rec("hRecEta_Rec4", nbinsEta_left, binsEta_left );
	TH1 *hRecPt_Rec4 = makeTH1Rec("hRecPt_Rec4", nbinsPt_left, binsPt_left );
	//For D0
	TH1 *hRecD0Eta_MC = makeTH1Rec("hRecD0Eta_MC", nbinsEta_left, binsEta_left );
	TH1 *hRecD0Pt_MC = makeTH1Rec("hRecD0Pt_MC", nbinsPt_left, binsPt_left );

	TH1 *hRecD0Eta_Rec1 = makeTH1Rec("hRecD0Eta_Rec1", nbinsEta_left, binsEta_left );
	TH1 *hRecD0Pt_Rec1 = makeTH1Rec("hRecD0Pt_Rec1", nbinsPt_left, binsPt_left );

	TH1 *hRecD0Eta_Rec2 = makeTH1Rec("hRecD0Eta_Rec2", nbinsEta_left, binsEta_left );
	TH1 *hRecD0Pt_Rec2 = makeTH1Rec("hRecD0Pt_Rec2", nbinsPt_left, binsPt_left );

	TH1 *hRecD0Eta_Rec3 = makeTH1Rec("hRecD0Eta_Rec3", nbinsEta_left, binsEta_left );
	TH1 *hRecD0Pt_Rec3 = makeTH1Rec("hRecD0Pt_Rec3", nbinsPt_left, binsPt_left );

	TH1 *hRecD0Eta_Rec4 = makeTH1Rec("hRecD0Eta_Rec4", nbinsEta_left, binsEta_left );
	TH1 *hRecD0Pt_Rec4 = makeTH1Rec("hRecD0Pt_Rec4", nbinsPt_left, binsPt_left );
#endif
	
	
	//for D0
	/*TH1 *hRecD0Eta_MC = makeTH1("hRecD0Eta_MC", 10,0,2.5, "hRecD0Eta_MC ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_MC = makeTH1("hRecD0Phi_MC", 10,0,5, "hRecD0Phi_MC ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_MC = makeTH1("hRecD0Pt_MC", 10,0,40, "hRecD0Pt_MC ; pT [GeV] ; Events ", kRed);

	TH1 *hRecD0Eta_Rec1 = makeTH1("hRecD0Eta_Rec1", 10,0,2.5, "hRecD0Eta_Rec1 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec1 = makeTH1("hRecD0Phi_Rec1", 10,0,5, "hRecD0Phi_Rec1 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec1 = makeTH1("hRecD0Pt_Rec1", 10,0,40, "hRecD0Pt_Rec1 ; pT [GeV] ; Events ", kRed);
	
	TH1 *hRecD0Eta_Rec2 = makeTH1("hRecD0Eta_Rec2", 10,0,2.5, "hRecD0Eta_Rec2 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec2 = makeTH1("hRecD0Phi_Rec2", 10,0,5, "hRecD0Phi_Rec2 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec2 = makeTH1("hRecD0Pt_Rec2", 10,0,40, "hRecD0Pt_Rec2 ; pT [GeV] ; Events ", kRed);
	
	TH1 *hRecD0Eta_Rec3 = makeTH1("hRecD0Eta_Rec3", 10,0,2.5, "hRecD0Eta_Rec3 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec3 = makeTH1("hRecD0Phi_Rec3", 10,0,5, "hRecD0Phi_Rec3 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec3 = makeTH1("hRecD0Pt_Rec3", 10,0,40, "hRecD0Pt_Rec3 ; pT [GeV] ; Events ", kRed);

	TH1 *hRecD0Eta_Rec4 = makeTH1("hRecD0Eta_Rec4", 10,0,2.5, "hRecD0Eta_Rec4 ; #eta ; Events ", kRed);
	TH1 *hRecD0Phi_Rec4 = makeTH1("hRecD0Phi_Rec4", 10,0,5, "hRecD0Phi_Rec4 ; #phi ; Events ", kRed);
	TH1 *hRecD0Pt_Rec4 = makeTH1("hRecD0Pt_Rec4", 10,0,40, "hRecD0Pt_Rec4 ; pT [GeV] ; Events ", kRed);*/

	//--------------------------------------------------
	//ADDRESSING THE MEMORY TO VECTOR AND VARIABLES
	//--------------------------------------------------
	TBranch *b_NameOfFiredTriggers; t1->SetBranchAddress("NameOfFiredTriggers",&NameOfFiredTriggers,&b_NameOfFiredTriggers);
#if GENvaraibles == 0	
#elif GENvaraibles == 1
	TBranch *b_PUWeight; t1->SetBranchAddress("PUWeight",&PUWeight,&b_PUWeight);
#endif
	//For D* 
	//-----------------------------------------
	//TBranch *b_D0mass = t1->GetBranch("D0mass");b_D0mass->SetAddress(&D0mass);
	TBranch *b_FlagDstarfromB; t1->SetBranchAddress("FlagDstarfromB",&FlagDstarfromB,&b_FlagDstarfromB);
	TBranch *b_Dsmass; t1->SetBranchAddress("Dsmass",&Dsmass,&b_Dsmass);
	TBranch *b_Dspt ; t1->SetBranchAddress("Dspt",&Dspt,&b_Dspt);
	TBranch *b_Dseta; t1->SetBranchAddress("Dseta",&Dseta,&b_Dseta);
	TBranch *b_Dsphi ; t1->SetBranchAddress("Dsphi",&Dsphi,&b_Dsphi);
	TBranch *b_DSDeltaR ; t1->SetBranchAddress("DSDeltaR",&DSDeltaR,&b_DSDeltaR);

	TBranch *b_D0mass; t1->SetBranchAddress("D0mass",&D0mass,&b_D0mass);
	TBranch *b_D0pt; t1->SetBranchAddress("D0pt",&D0pt,&b_D0pt);
	TBranch *b_D0eta; t1->SetBranchAddress("D0eta",&D0eta,&b_D0eta);
	TBranch *b_D0phi; t1->SetBranchAddress("D0phi",&D0phi,&b_D0phi);
	TBranch *b_Anglephi ; t1->SetBranchAddress("Anglephi",&Anglephi,&b_Anglephi);
	TBranch *b_D0fromDSsXY ; t1->SetBranchAddress("D0fromDSsXY",&D0fromDSsXY,&b_D0fromDSsXY);
	TBranch *b_D0fromDSdXY ; t1->SetBranchAddress("D0fromDSdXY",&D0fromDSdXY,&b_D0fromDSdXY);
	TBranch *b_D0fromDSs3D ; t1->SetBranchAddress("D0fromDSs3D",&D0fromDSs3D,&b_D0fromDSs3D);
	TBranch *b_D0fromDSd3D ; t1->SetBranchAddress("D0fromDSd3D",&D0fromDSd3D,&b_D0fromDSd3D);
	TBranch *b_D0_VtxProb ; t1->SetBranchAddress("D0_VtxProb",&D0_VtxProb,&b_D0_VtxProb);

	/*TBranch *b_TrkSpt ; t1->SetBranchAddress("TrkSpt",&TrkSpt,&b_TrkSpt);
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
	TBranch *b_TrkKphi ; t1->SetBranchAddress("TrkKphi",&TrkKphi,&b_TrkKphi);*/
	//----------------------------------------
	//For D* WRONG COMBINATION
	//-----------------------------------------
	TBranch *b_DsmassWrong; t1->SetBranchAddress("DsmassWrong",&DsmassWrong,&b_DsmassWrong);
	TBranch *b_DsetaWrong ; t1->SetBranchAddress("DsetaWrong",&DsetaWrong,&b_DsetaWrong);
	TBranch *b_DsphiWrong ; t1->SetBranchAddress("DsphiWrong",&DsphiWrong,&b_DsphiWrong);
	TBranch *b_DsptWrong ; t1->SetBranchAddress("DsptWrong",&DsptWrong,&b_DsptWrong);

	TBranch *b_D0massWrong; t1->SetBranchAddress("D0massWrong",&D0massWrong,&b_D0massWrong);
	TBranch *b_D0etaWrong; t1->SetBranchAddress("D0etaWrong",&D0etaWrong,&b_D0etaWrong);
	TBranch *b_D0phiWrong; t1->SetBranchAddress("D0phiWrong",&D0phiWrong,&b_D0phiWrong);
	TBranch *b_D0ptWrong; t1->SetBranchAddress("D0ptWrong",&D0ptWrong,&b_D0ptWrong);
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
	TBranch *b_D0KpidXY ; t1->SetBranchAddress("D0KpidXY",&D0KpidXY,&b_D0KpidXY);
	TBranch *b_D0KpisXY ; t1->SetBranchAddress("D0KpisXY",&D0KpisXY,&b_D0KpisXY);
	TBranch *b_D0Kpid3D ; t1->SetBranchAddress("D0Kpid3D",&D0Kpid3D,&b_D0Kpid3D);
	TBranch *b_D0Kpie3D ; t1->SetBranchAddress("D0Kpie3D",&D0Kpie3D,&b_D0Kpie3D);
	TBranch *b_D0_kT ; t1->SetBranchAddress("D0KpikT",&D0_kT,&b_D0_kT);

	/*TBranch *b_TrkD0Kdxy ; t1->SetBranchAddress("TrkD0Kdxy",&TrkD0Kdxy,&b_TrkD0Kdxy);
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
	TBranch *b_TrkD0pinhits ; t1->SetBranchAddress("TrkD0pinhits",&TrkD0pinhits,&b_TrkD0pinhits);*/

#if GENvaraibles == 0	
#elif GENvaraibles == 1
	//----------------------------------------
	//FOR D* MC
	//-----------------------------------------
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
	TBranch *b_MCD0lifetime ; t1->SetBranchAddress("MCD0lifetime",&MCD0lifetime,&b_MCD0lifetime);

	/*TBranch *b_MCDsKphi ; t1->SetBranchAddress("MCDsKphi",&MCDsKphi,&b_MCDsKphi);
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
	TBranch *b_MCDsPimass ; t1->SetBranchAddress("MCDsPimass",&MCDsPimass,&b_MCDsPimass);*/
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
	TBranch *b_MCpromptD0lifetime ; t1->SetBranchAddress("MCpromptD0lifetime",&MCpromptD0lifetime,&b_MCpromptD0lifetime); 
	TBranch *b_MCpromptD0_DispAngle ; t1->SetBranchAddress("MCpromptD0_DispAngle",&MCpromptD0_DispAngle,&b_MCpromptD0_DispAngle);

	/*TBranch *b_MCpromptD0_Keta ; t1->SetBranchAddress("MCpromptD0_Keta",&MCpromptD0_Keta,&b_MCpromptD0_Keta);
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
	TBranch *b_MCpromptD0_Pimass ; t1->SetBranchAddress("MCpromptD0_Pimass",&MCpromptD0_Pimass,&b_MCpromptD0_Pimass);*/
#endif
	//TBranch *b_teste = t1->GetBranch("teste");b_teste->SetAddress(&teste);	

	for (Long64_t jentry=0; jentry < nentries; jentry++) //loop tree entries for file f1
	{
		Long64_t ientry = t1->LoadTree(jentry);
      //cout << "nentries " << nentries << " ientry " << ientry << " jentry " << jentry <<endl;
      if (ientry < 0) break;

		NumberOfEvents++;

		//Output about percent program executed
		double percent = (jentry*100)/nentries;		
		if ( jentry % partpercent == 0) cout <<"*****"<< percent << "per cent done***" << endl;
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
		//	cout << "NameOfFiredTriggers ["<< tescont <<"]: " << teststring  << endl;}
#if GENvaraibles == 0	
#elif GENvaraibles == 1
		b_PUWeight->GetEntry(ientry);
		//cout << "PUWeight: " << PUWeight << endl;
#endif

		EventsAfterTrigger++;
		//----------------------------------------
		//For D* 
		//-----------------------------------------
		if ( FlagDs){		
		b_FlagDstarfromB->GetEntry(ientry);
		b_Dsmass->GetEntry(ientry);
		b_Dseta->GetEntry(ientry);
		b_Dsphi->GetEntry(ientry);
		b_Dspt->GetEntry(ientry);

		b_D0mass->GetEntry(ientry);
		b_D0eta->GetEntry(ientry);
		b_D0phi->GetEntry(ientry);
		b_D0pt->GetEntry(ientry);
		b_Anglephi->GetEntry(ientry);
		b_D0fromDSsXY->GetEntry(ientry);
		b_D0fromDSdXY->GetEntry(ientry);
		b_D0fromDSs3D->GetEntry(ientry); 
		b_D0fromDSd3D->GetEntry(ientry);
		b_D0_VtxProb->GetEntry(ientry);
		b_DSDeltaR->GetEntry(ientry);

		/*b_TrkSpt->GetEntry(ientry);
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
		b_Trkpiphi->GetEntry(ientry);*/
		}
		//----------------------------------------
		//For D* Wrong combination
		//-----------------------------------------
		if( FlagDsWrong ){
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
		}
		//----------------------------------------
		//For D0 
		//-----------------------------------------
		if ( FlagD0){
		b_D0Kpimass->GetEntry(ientry);
		b_D0Kpi_VtxProb->GetEntry(ientry);
		b_D0Kpis3D->GetEntry(ientry);
		b_D0KpiDispAngle->GetEntry(ientry);
		b_D0KpidXY->GetEntry(ientry);
		b_D0KpisXY->GetEntry(ientry);
		b_D0Kpid3D->GetEntry(ientry);
		b_D0Kpie3D->GetEntry(ientry);
		b_D0_kT->GetEntry(ientry);

		b_D0Kpipt->GetEntry(ientry);
		b_D0Kpieta->GetEntry(ientry);
		b_D0Kpiphie->GetEntry(ientry);
	
		/*b_TrkD0Kdxy->GetEntry(ientry);
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
		b_TrkD0pinhits->GetEntry(ientry);*/
		}
		//----------------------------------------
		//For D* MC
		//-----------------------------------------
#if GENvaraibles == 0
#elif GENvaraibles == 1
		if( FlagDsMC){
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
		b_MCD0lifetime->GetEntry(ientry);


		/*b_MCDsKphi->GetEntry(ientry);
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
		b_MCDsPimass->GetEntry(ientry);*/
		}
		//----------------------------------------
		//For D0 MC
		//-----------------------------------------
		if( FlagD0MC){
		b_MCpromptD0eta->GetEntry(ientry);
		b_MCpromptD0phi->GetEntry(ientry);
		b_MCpromptD0pt->GetEntry(ientry);
		b_MCpromptD0energy->GetEntry(ientry);
		b_MCpromptD0p->GetEntry(ientry);
		b_MCpromptD0et->GetEntry(ientry);
		b_MCpromptD0rapidity->GetEntry(ientry);
		b_MCpromptD0mass->GetEntry(ientry); 
		b_MCpromptD0_DispAngle->GetEntry(ientry);
		b_MCpromptD0lifetime->GetEntry(ientry);

		/*b_MCpromptD0_Keta->GetEntry(ientry);
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
		b_MCpromptD0_Pimass->GetEntry(ientry);*/
		}
#endif
		//----------------------------------------
		//For D* and D0(from D*)
		//-----------------------------------------
		if ( (D0mass->size() > 0) && FlagDs){
		for(unsigned int i=0; i < D0mass->size(); i++)
		{  
			PossibleDs++;					
			if( D0_VtxProb->at(i) < 0.01) continue; DsD0ProbVtx++;		
			if( Anglephi->at(i) < 0.99 ) continue; DsD0Angle++;	
			if( fabs(D0mass->at(i) - 1.86484) > 0.1 ) continue; D0fromDsMinusPDG++;
			if ( D0pt->at(i) < 3. ) continue; DsD0pt++;
			if( (Dsmass->at(i) - D0mass->at(i)) > 0.16) continue; DsD0deltaM++;
			DsmassHisto0->Fill(Dsmass->at(i));
			TotalD0fromDSs3DHisto->Fill(D0fromDSs3D->at(i));
			if( D0fromDSs3D->at(i)>0.5) { DsmassHisto1->Fill(Dsmass->at(i));}
			if( D0fromDSs3D->at(i)>1.0) { DsmassHisto2->Fill(Dsmass->at(i));}
			if( D0fromDSs3D->at(i)>1.5) { DsmassHisto3->Fill(Dsmass->at(i));}
			if( D0fromDSs3D->at(i)>2.0) { DsmassHisto4->Fill(Dsmass->at(i));}
			if( D0fromDSs3D->at(i)>2.5) { DsmassHisto5->Fill(Dsmass->at(i));}
			if( D0fromDSs3D->at(i)>3.0) { DsmassHisto6->Fill(Dsmass->at(i));}

			if( D0fromDSs3D->at(i) < SigCutD0fromDs) continue; DsD0Sig++;

			if( D0pt->at(i)>0.0) { DsmassPtHisto1->Fill(Dsmass->at(i));}
			if( D0pt->at(i)>0.5) { DsmassPtHisto2->Fill(Dsmass->at(i));}
			if( D0pt->at(i)>1.0) { DsmassPtHisto3->Fill(Dsmass->at(i));}
			if( D0pt->at(i)>1.5) { DsmassPtHisto4->Fill(Dsmass->at(i));}
			if( D0pt->at(i)>2.0) { DsmassPtHisto5->Fill(Dsmass->at(i));}
			if( D0pt->at(i)>2.5) { DsmassPtHisto6->Fill(Dsmass->at(i));}
			if( D0pt->at(i)>3.0) { DsmassPtHisto7->Fill(Dsmass->at(i));}
			if( D0pt->at(i)>3.5) { DsmassPtHisto8->Fill(Dsmass->at(i));}
			if( D0pt->at(i)>4.0) { DsmassPtHisto9->Fill(Dsmass->at(i));}		
		
			D0massHisto->Fill(D0mass->at(i));
			DsmassHisto->Fill(Dsmass->at(i));
			DsMinusD0Histo->Fill((Dsmass->at(i) - D0mass->at(i)));
			D0fromDSs3DHisto->Fill(D0fromDSs3D->at(i));
			D0VtxProbHisto->Fill(D0_VtxProb->at(i));
			DslifetimeHisto->Fill( (D0mass->at(i)*D0fromDSdXY->at(i)*pow(10,-2)) / (c*D0pt->at(i)) );

			var_Dsmass = Dsmass->at(i);
			var_Dseta = Dseta->at(i);
			var_Dsphi = Dsphi->at(i);
			var_Dspt = Dspt->at(i);
			var_Dslifetime = ( D0mass->at(i) * D0fromDSdXY->at(i)*pow(10,-2) ) / ( c * D0pt->at(i) );
			var_DSDeltaR = DSDeltaR->at(i);

			var_DsMinusD0 = (Dsmass->at(i) - D0mass->at(i));
			var_D0mass = D0mass->at(i);
			var_D0eta = D0eta->at(i);
			var_D0phi = D0phi->at(i);
			var_D0pt = D0pt->at(i);
			var_D0fromDSsXY = D0fromDSsXY->at(i);
			var_D0fromDSs3D = D0fromDSs3D->at(i);
			var_Anglephi = Anglephi->at(i);
			var_D0_VtxProb = D0_VtxProb->at(i);

			/*var_TrkKpt = TrkKpt->at(i);		
			var_TrkKnhits = TrkKnhits->at(i);
			var_TrkKchi2 = TrkKchi2->at(i);
			var_TrkKdxy = TrkKdxy->at(i);
			var_TrkKdz = TrkKdz->at(i);
			var_TrkKeta = TrkKeta->at(i);
			var_TrkKphi = TrkKphi->at(i);

			var_Trkpipt = Trkpipt->at(i);
			var_Trkpinhits = Trkpinhits->at(i);
			var_Trkpichi2 = Trkpichi2->at(i);
			var_Trkpidxy = Trkpidxy->at(i);
			var_Trkpidz = Trkpidz->at(i);
			var_Trkpieta = Trkpieta->at(i);
			var_Trkpiphi = Trkpiphi->at(i);

			var_TrkSpt = TrkSpt->at(i);
			var_TrkSchi2 = TrkSchi2->at(i);
			var_TrkSnhits = TrkSnhits->at(i);
			var_TrkSdz = TrkSdz->at(i);
			var_TrkSdxy = TrkSdxy->at(i);
			var_TrkSeta = TrkSeta->at(i);
			var_TrkSphi = TrkSphi->at(i);*/
			//var_TrkScharge = TrkScharge->at(i);
				
#if GENvaraibles == 0
			//Histograms
			hDsmass->Fill(Dsmass->at(i));
			hDseta->Fill(Dseta->at(i));
			hDsphi->Fill(Dsphi->at(i));
			hDspt->Fill(Dspt->at(i));
			hDslifetime->Fill(( D0mass->at(i)*D0fromDSdXY->at(i)*pow(10,-2))/(c*D0pt->at(i)));
			hDsMinusD0->Fill(( Dsmass->at(i) - D0mass->at(i)));

			hD0mass->Fill(D0mass->at(i));
			hD0eta->Fill(D0eta->at(i));
			hD0phi->Fill(D0phi->at(i));
			hD0pt->Fill(D0pt->at(i));

			/*hTrkKpt->Fill(TrkKpt->at(i));
			hTrkKchi2->Fill(TrkKchi2->at(i));
			hTrkKnhits->Fill(TrkKnhits->at(i));
			hTrkKdz->Fill(TrkKdz->at(i));
			hTrkKdxy->Fill(TrkKdxy->at(i));
			hTrkKeta->Fill(TrkKeta->at(i));
			hTrkKphi->Fill(TrkKphi->at(i));

			hTrkpipt->Fill(Trkpipt->at(i));
			hTrkpichi2->Fill(Trkpichi2->at(i));
			hTrkpinhits->Fill(Trkpinhits->at(i));
			hTrkpidz->Fill(Trkpidz->at(i));
			hTrkpidxy->Fill(Trkpidxy->at(i));
			hTrkpieta->Fill(Trkpieta->at(i));
			hTrkpiphi->Fill(Trkpiphi->at(i));

			hTrkSpt->Fill(TrkSpt->at(i));
			hTrkSchi2->Fill(TrkSchi2->at(i));
			hTrkSnhits->Fill(TrkSnhits->at(i));
			hTrkSdz->Fill(TrkSdz->at(i));
			hTrkSdxy->Fill(TrkSdxy->at(i));
			hTrkSeta->Fill(TrkSeta->at(i));
			hTrkSphi->Fill(TrkSphi->at(i));*/
#elif GENvaraibles == 1
			//Histograms
			hDsmass->Fill(Dsmass->at(i),PUWeight);
			hDseta->Fill(Dseta->at(i),PUWeight);
			hDsphi->Fill(Dsphi->at(i),PUWeight);
			hDspt->Fill(Dspt->at(i),PUWeight);
			hDslifetime->Fill((D0mass->at(i)*D0fromDSdXY->at(i)*pow(10,-2))/(c*D0pt->at(i)),PUWeight);
			hDsMinusD0->Fill((Dsmass->at(i) - D0mass->at(i)),PUWeight);

			hD0mass->Fill(D0mass->at(i),PUWeight);
			hD0eta->Fill(D0eta->at(i),PUWeight);
			hD0phi->Fill(D0phi->at(i),PUWeight);
			hD0pt->Fill(D0pt->at(i),PUWeight);

			/*hTrkKpt->Fill(TrkKpt->at(i),PUWeight);
			hTrkKchi2->Fill(TrkKchi2->at(i),PUWeight);
			hTrkKnhits->Fill(TrkKnhits->at(i),PUWeight);
			hTrkKdz->Fill(TrkKdz->at(i),PUWeight);
			hTrkKdxy->Fill(TrkKdxy->at(i),PUWeight);
			hTrkKeta->Fill(TrkKeta->at(i),PUWeight);
			hTrkKphi->Fill(TrkKphi->at(i),PUWeight);

			hTrkpipt->Fill(Trkpipt->at(i),PUWeight);
			hTrkpichi2->Fill(Trkpichi2->at(i),PUWeight);
			hTrkpinhits->Fill(Trkpinhits->at(i),PUWeight);
			hTrkpidz->Fill(Trkpidz->at(i),PUWeight);
			hTrkpidxy->Fill(Trkpidxy->at(i),PUWeight);
			hTrkpieta->Fill(Trkpieta->at(i),PUWeight);
			hTrkpiphi->Fill(Trkpiphi->at(i),PUWeight);

			hTrkSpt->Fill(TrkSpt->at(i),PUWeight);
			hTrkSchi2->Fill(TrkSchi2->at(i),PUWeight);
			hTrkSnhits->Fill(TrkSnhits->at(i),PUWeight);
			hTrkSdz->Fill(TrkSdz->at(i),PUWeight);
			hTrkSdxy->Fill(TrkSdxy->at(i),PUWeight);
			hTrkSeta->Fill(TrkSeta->at(i),PUWeight);
			hTrkSphi->Fill(TrkSphi->at(i),PUWeight);*/

			//----------------------------------------
			//For D* Contamination
			//-----------------------------------------
			DstarTotal++;
			if ( FlagDstarfromB->size() != 0){ //Protection
				for(unsigned int j = 0; j < FlagDstarfromB->size(); j++)
				{  													
					if( FlagDstarfromB->at(j) == 1) { DstarSec++;}
				}
			} //Protection
#endif
			
			t_analysis.Fill();

			if (debug)cout << "debug D* 12 --------------------"<< endl;
  		}
		}
		//----------------------------------------
		//For D* Wrong Combination
		//-----------------------------------------
		if( (D0massWrong->size() > 0) && FlagDsWrong){
		for(unsigned int i=0; i < D0massWrong->size(); i++)
		{ 	if (debug)cout << "debug D* 13 --------------------"<< endl; 			
			if( D0_VtxProbWrong->at(i) < 0.01) continue;
			if( AnglephiWrong->at(i) < 0.99 ) continue;
			if( fabs(D0massWrong->at(i) - 1.86484) > 0.1 ) continue;
			if ( D0ptWrong->at(i) < 3. ) continue;
			if( (DsmassWrong->at(i) - D0massWrong->at(i)) > 0.16) continue;
			if( D0fromDSs3DWrong->at(i) < SigCutD0fromDs) continue;
					
			var_DsMinusD0 = (DsmassWrong->at(i) - D0massWrong->at(i));
			var_D0massWrong = D0massWrong->at(i);
			var_DsmassWrong = DsmassWrong->at(i);
			var_D0etaWrong = D0etaWrong->at(i);
			var_D0phiWrong = D0phiWrong->at(i);
			var_DsetaWrong = DsetaWrong->at(i);
			var_DsphiWrong = DsphiWrong->at(i);
			var_D0ptWrong = D0ptWrong->at(i);
			var_DsptWrong = DsptWrong->at(i);

			t_DstarWrongCombination.Fill();
			if (debug)cout << "debug D* 14 --------------------"<< endl;
  		}//For D* Wrong Combination
		}
		//----------------------------------------
		//For D0 
		//-----------------------------------------
		if (D0Kpimass->size() > 0 && FlagD0){
		for(unsigned int i=0; i < D0Kpimass->size(); i++)
		{  if (debug)cout << "debug D* 15 --------------------"<< endl;
			PossibleD0++;
			if( D0Kpi_VtxProb->at(i) < 0.01) continue; D0ProbVtx++;	
			if( D0KpiDispAngle->at(i) < 0.99 ) continue; D0Angle++;
			TotalD0Kpis3DHistoHisto->Fill(D0Kpis3D->at(i));
			if( abs(D0Kpimass->at(i)-1.86484) > 0.1) continue;	CountD0minusPDG++;
			
			if( D0Kpis3D->at(i) < SigCutD0)
			{	//For study a pT cut in the mass
				if( D0Kpipt->at(i)>0.0) { D0KpimassPtHisto1->Fill(D0Kpimass->at(i));}
				if( D0Kpipt->at(i)>0.5) { D0KpimassPtHisto2->Fill(D0Kpimass->at(i));}
				if( D0Kpipt->at(i)>1.0) { D0KpimassPtHisto3->Fill(D0Kpimass->at(i));}
				if( D0Kpipt->at(i)>1.5) { D0KpimassPtHisto4->Fill(D0Kpimass->at(i));}
				if( D0Kpipt->at(i)>2.0) { D0KpimassPtHisto5->Fill(D0Kpimass->at(i));}
				if( D0Kpipt->at(i)>2.5) { D0KpimassPtHisto6->Fill(D0Kpimass->at(i));}
				if( D0Kpipt->at(i)>3.0) { D0KpimassPtHisto7->Fill(D0Kpimass->at(i));}
				if( D0Kpipt->at(i)>3.5) { D0KpimassPtHisto8->Fill(D0Kpimass->at(i));}
				if( D0Kpipt->at(i)>4.0) { D0KpimassPtHisto9->Fill(D0Kpimass->at(i));}
			}

			if( D0Kpipt->at(i) < 2. ) continue; CountD0pt++;

			//For study a Significance cut in the mass
			D0KpimassHisto0->Fill(D0Kpimass->at(i));
			if( D0Kpis3D->at(i)>0.5) { D0KpimassHisto1->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>1.0) { D0KpimassHisto2->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>1.5) { D0KpimassHisto3->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>2.0) { D0KpimassHisto4->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>2.5) { D0KpimassHisto5->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>3.0) { D0KpimassHisto6->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>3.5) { D0KpimassHisto7->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>4.0) { D0KpimassHisto8->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>4.5) { D0KpimassHisto9->Fill(D0Kpimass->at(i));}
			if( D0Kpis3D->at(i)>5.0) { D0KpimassHisto10->Fill(D0Kpimass->at(i));}

			if( D0Kpis3D->at(i) < SigCutD0) continue; D0Sig++;
			
			//Fill Variables
			D0KpimassHisto->Fill(D0Kpimass->at(i));
			D0Kpis3DHisto->Fill(D0Kpis3D->at(i));
			D0Kpi_VtxProbHisto->Fill(D0Kpi_VtxProb->at(i));
			D0lifetimeHisto->Fill((D0Kpimass->at(i)*D0KpidXY->at(i)*pow(10,-2))/(c*D0Kpipt->at(i)));

			var_D0Kpimass = D0Kpimass->at(i);
			var_D0Kpi_VtxProb = D0Kpi_VtxProb->at(i);
			var_D0Kpis3D = D0Kpis3D->at(i);
			var_D0KpiDispAngle = D0KpiDispAngle->at(i);
			var_D0KpisXY = D0KpisXY->at(i);
			var_D0KpidXY = D0KpidXY->at(i);
			var_D0Kpie3D = D0Kpie3D->at(i);
			var_D0_kT = D0_kT->at(i);
			var_D0Kpipt = D0Kpipt->at(i);
			var_D0Kpieta = D0Kpieta->at(i);
			var_D0Kpis3D = D0Kpis3D->at(i);
			var_D0Kpiphi = D0Kpiphi->at(i);
			var_D0lifetime = (D0Kpimass->at(i)*D0KpidXY->at(i)*pow(10,-2))/(c*D0Kpipt->at(i));

			/*var_TrkD0Kdxy = TrkD0Kdxy->at(i);
			var_TrkD0Kdz = TrkD0Kdz->at(i);
			var_TrkD0Kchi2 = TrkD0Kchi2->at(i);
			var_TrkD0Kpt = TrkD0Kpt->at(i);
			var_TrkD0Keta = TrkD0Keta->at(i);
			var_TrkD0kphi = TrkD0kphi->at(i);
			var_TrkD0Knhits = TrkD0Knhits->at(i);

			var_TrkD0pidxy = TrkD0pidxy->at(i);
			var_TrkD0pidz = TrkD0pidz->at(i);
			var_TrkD0pichi2 = TrkD0pichi2->at(i);
			var_TrkD0pieta = TrkD0pieta->at(i);
			var_TrkD0piphi = TrkD0piphi->at(i);
			var_TrkD0pinhits = TrkD0pinhits->at(i);*/
#if GENvaraibles == 0
			//Fill histograms
			hD0Kpimass->Fill(D0Kpimass->at(i));
			hD0Kpieta->Fill(D0Kpieta->at(i));
			hD0Kpiphi->Fill(D0Kpiphi->at(i));
			hD0Kpipt->Fill(D0Kpipt->at(i));
			hD0lifetime->Fill((D0Kpimass->at(i)*D0KpidXY->at(i)*pow(10,-2))/(c*D0Kpipt->at(i)));

			/*hTrkD0Kdxy->Fill(TrkD0Kdxy->at(i));
			hTrkD0Kdz->Fill(TrkD0Kdz->at(i));
			hTrkD0Kchi2->Fill(TrkD0Kchi2->at(i));
			hTrkD0Kpt->Fill(TrkD0Kpt->at(i));
			hTrkD0Keta->Fill(TrkD0Keta->at(i));
			hTrkD0kphi->Fill(TrkD0kphi->at(i));
			hTrkD0Knhits->Fill(TrkD0Knhits->at(i));

			hTrkD0pidxy->Fill(TrkD0pidxy->at(i));
			hTrkD0pidz->Fill(TrkD0pidz->at(i));
			hTrkD0pichi2->Fill(TrkD0pichi2->at(i));
			hTrkD0pipt->Fill(TrkD0pipt->at(i));
			hTrkD0pieta->Fill(TrkD0pieta->at(i));
			hTrkD0piphi->Fill(TrkD0piphi->at(i));
			hTrkD0pinhits->Fill(TrkD0pinhits->at(i));*/
#elif GENvaraibles == 1
			//histogramas
			hD0Kpimass->Fill(D0Kpimass->at(i),PUWeight);
			hD0Kpieta->Fill(D0Kpieta->at(i),PUWeight);
			hD0Kpiphi->Fill(D0Kpiphi->at(i),PUWeight);
			hD0Kpipt->Fill(D0Kpipt->at(i),PUWeight);
			hD0lifetime->Fill((D0Kpimass->at(i)*D0KpidXY->at(i)*pow(10,-2))/(c*D0Kpipt->at(i)),PUWeight);

			/*hTrkD0Kdxy->Fill(TrkD0Kdxy->at(i),PUWeight);
			hTrkD0Kdz->Fill(TrkD0Kdz->at(i),PUWeight);
			hTrkD0Kchi2->Fill(TrkD0Kchi2->at(i),PUWeight);
			hTrkD0Kpt->Fill(TrkD0Kpt->at(i),PUWeight);
			hTrkD0Keta->Fill(TrkD0Keta->at(i),PUWeight);
			hTrkD0kphi->Fill(TrkD0kphi->at(i),PUWeight);
			hTrkD0Knhits->Fill(TrkD0Knhits->at(i),PUWeight);

			hTrkD0pidxy->Fill(TrkD0pidxy->at(i),PUWeight);
			hTrkD0pidz->Fill(TrkD0pidz->at(i),PUWeight);
			hTrkD0pichi2->Fill(TrkD0pichi2->at(i),PUWeight);
			hTrkD0pipt->Fill(TrkD0pipt->at(i),PUWeight);
			hTrkD0pieta->Fill(TrkD0pieta->at(i),PUWeight);
			hTrkD0piphi->Fill(TrkD0piphi->at(i),PUWeight);
			hTrkD0pinhits->Fill(TrkD0pinhits->at(i),PUWeight);*/
#endif
			//vectorInvariantMass_D0.push_back(D0mass->at(i));

			t_D0analysis.Fill();	
  		}//For D0(Direct)
		}

#if GENvaraibles == 0	
#elif GENvaraibles == 1
		//----------------------------------------
		//For D* MC
		//-----------------------------------------	
		//vector<double> deltaR_vec;
		if ( (MCDsmass->size() > 0) && FlagDsMC){ //D* MC protection
		for(unsigned int i=0; i < MCDsmass->size(); i++)
		{  if (debug)cout << "For D* MC loop --------------------"<< endl;
			if ( MCDsmass->size() == 0) continue;	
			/*var_MCDseta = MCDseta->at(i);
			var_MCDsphi = MCDsphi->at(i);
			var_MCDspt = MCDspt->at(i);
			var_MCDsenergy = MCDsenergy->at(i);
			var_MCDsp = MCDsp->at(i);
			var_MCDset = MCDset->at(i);
			var_MCDsmass = MCDsmass->at(i);

			var_MCD0eta = MCD0eta->at(i);
			var_MCD0phi = MCD0phi->at(i);
			var_MCD0pt = MCD0pt->at(i);
			var_MCD0energy = MCD0energy->at(i);
			var_MCD0p = MCD0p->at(i);
			var_MCD0et = MCD0et->at(i);
			var_MCD0rapidity = MCD0rapidity->at(i);
			var_MCD0mass = MCD0mass->at(i);*/
			var_MCD0lifetime = MCD0lifetime->at(i);

			/*var_MCDsKphi = MCDsKphi->at(i);
			var_MCDsKpt = MCDsKpt->at(i);
			var_MCDsKenergy = MCDsKenergy->at(i);
			var_MCDsKp = MCDsKp->at(i);
			var_MCDsKet = MCDsKet->at(i);
			var_MCDsKrapidity = MCDsKrapidity->at(i);
			var_MCDsKmass = MCDsKmass->at(i);

			var_MCDsPieta = MCDsPieta->at(i);
			var_MCDsPiphi = MCDsPiphi->at(i);
			var_MCDsPipt = MCDsPipt->at(i);
			var_MCDsPienergy = MCDsPienergy->at(i);
			var_MCDsPip = MCDsPip->at(i);
			var_MCDsPiet = MCDsPiet->at(i);
			var_MCDsPirapidity = MCDsPirapidity->at(i);
			var_MCDsPimass = MCDsPimass->at(i);*/

			hRecEta_MC->Fill(abs(MCDseta->at(i)));
			hRecPt_MC->Fill(MCDspt->at(i));

			//for(unsigned int i=0; i < MCDsmass->size(); i++)
			 

			if( (MCDspt->at(i) > 4.) and (MCDspt->at(i) < 5.) ){ GenDspT1++;}
			if( (MCDspt->at(i) > 5.) and (MCDspt->at(i) < 6.) ){ GenDspT2++;}
			if( (MCDspt->at(i) > 6.) and (MCDspt->at(i) < 7.) ){ GenDspT3++;}
			if( (MCDspt->at(i) > 7.) and (MCDspt->at(i) < 8.) ){ GenDspT4++;}
			if( (MCDspt->at(i) > 8.) and (MCDspt->at(i) < 12.) ){ GenDspT5++;}
			if( (MCDspt->at(i) > 12.) and (MCDspt->at(i) < 16.) ){ GenDspT6++;}
			if( (MCDspt->at(i) > 16.) and (MCDspt->at(i) < 24.) ){ GenDspT7++;}
			if( (MCDspt->at(i) > 24.) and (MCDspt->at(i) < 40.) ){ GenDspT8++;}
			if( (MCDspt->at(i) > 40.) and (MCDspt->at(i) < 100.) ){ GenDspT9++;}

			if( (abs(MCDseta->at(i)) > 0.) and (abs(MCDseta->at(i)) < 0.2) ){ GenDsEta1++;}
			if( (abs(MCDseta->at(i)) > 0.2) and (abs(MCDseta->at(i)) < 0.4) ){ GenDsEta2++;}
			if( (abs(MCDseta->at(i)) > 0.4) and (abs(MCDseta->at(i)) < 0.6) ){ GenDsEta3++;}
			if( (abs(MCDseta->at(i)) > 0.6) and (abs(MCDseta->at(i)) < 0.8) ){ GenDsEta4++;}
			if( (abs(MCDseta->at(i)) > 0.8) and (abs(MCDseta->at(i)) < 1.0) ){ GenDsEta5++;}
			if( (abs(MCDseta->at(i)) > 1.0) and (abs(MCDseta->at(i)) < 1.2) ){ GenDsEta6++;}
			if( (abs(MCDseta->at(i)) > 1.2) and (abs(MCDseta->at(i)) < 1.4) ){ GenDsEta7++;}
			if( (abs(MCDseta->at(i)) > 1.4) and (abs(MCDseta->at(i)) < 1.6) ){ GenDsEta8++;}
			if( (abs(MCDseta->at(i)) > 1.6) and (abs(MCDseta->at(i)) < 1.8) ){ GenDsEta9++;}
			if( (abs(MCDseta->at(i)) > 1.8) and (abs(MCDseta->at(i)) < 2.1) ){ GenDsEta10++;}

			//t_DsMC.Fill();
			if (debug)cout << "debug D* 18 --------------------"<< endl;
			//----------------------------------------
			//For Ds Matching
			//-----------------------------------------
			if ( D0mass->size() == 0) continue;
			double minDeltaR = 1000.;
			int id = -1;
			double deltaR = 0.;
			for(unsigned int j=0; j < D0mass->size(); j++)
			{  if (debug)cout << "debug D* 19 --------------------"<< endl;													
				if( D0_VtxProb->at(j) < 0.01) continue;
				if( Anglephi->at(j) < 0.99 ) continue;
				if( D0pt->at(j) < 3. ) continue;
				if( (Dsmass->at(j) - D0mass->at(j)) > 0.16) continue; 
				if( D0fromDSs3D->at(j) < SigCutD0fromDs) continue;
				//DeltaR Function
				deltaR = sqrt( pow( (MCDseta->at(i)-Dseta->at(j)) ,2) + pow( (MCDsphi->at(i)-Dsphi->at(j)) ,2));
				if( deltaR < minDeltaR) { id = j; minDeltaR = deltaR;} //Get the lower delta
				if (debug)cout << "debug D* 20 --------------------"<< endl;
			}	

			if (id == -1) continue; // protection
			var_MCD0lifetimeMatching = MCD0lifetime->at(i);
			/*var_MCDsetaMatching = MCDseta->at(i);
			var_MCDsphiMatching = MCDsphi->at(i);
			var_MCDsptMatching = MCDspt->at(i);

			var_DsptMatching = Dspt->at(id);
			var_DsetaMatching = Dseta->at(id);
			var_DsphiMatching = Dsphi->at(id);*/

			DsdeltaEta = MCDseta->at(i) - Dseta->at(id);
			DsdeltaPhi = MCDsphi->at(i) - Dsphi->at(id);
			DsdeltaPt = MCDspt->at(i) - Dspt->at(id);
			var_deltaRDsMatching = deltaR;

			double recMass = -1; double GenMass = -1; // Protection
			double recPt = -1; double GenPt = -1; // Protection
			double recEta = -1; double GenEta = -1; // Protection
			double recTime = -1; double GenTime = -1; // Protection
			//Histograms for efficiency
			if( deltaR < 0.03){
				hRecEta_Rec1->Fill(abs(Dseta->at(id)));
				hRecPt_Rec1->Fill(Dspt->at(id));

				recMass = Dsmass->at(id); GenMass = MCDsmass->at(i);
				recPt = Dspt->at(id); GenPt = MCDspt->at(i);
				recEta = Dseta->at(id); GenEta = MCDseta->at(i);
				recTime = (D0mass->at(id)*D0fromDSdXY->at(id)*pow(10,-2)) / (c*D0pt->at(id)); 
				GenTime = MCD0lifetime->at(i);

				hdeltaRDstar->Fill( recMass, GenMass );
				hdeltaRDstarPt->Fill( recPt, GenPt );
				hdeltaRDstarEta->Fill( abs(recEta), abs(GenEta) );
				hdeltaRDstarTime->Fill( recTime, GenTime );

				if( (Dspt->at(id) > 4.) and (Dspt->at(id) < 5.) ){ RecDspT1++;}
				if( (Dspt->at(id) > 5.) and (Dspt->at(id) < 6.) ){ RecDspT2++;}
				if( (Dspt->at(id) > 6.) and (Dspt->at(id) < 7.) ){ RecDspT3++;}
				if( (Dspt->at(id) > 7.) and (Dspt->at(id) < 8.) ){ RecDspT4++;}
				if( (Dspt->at(id) > 8.) and (Dspt->at(id) < 12.) ){ RecDspT5++;}
				if( (Dspt->at(id) > 12.) and (Dspt->at(id) < 16.) ){ RecDspT6++;}
				if( (Dspt->at(id) > 16.) and (Dspt->at(id) < 24.) ){ RecDspT7++;}
				if( (Dspt->at(id) > 24.) and (Dspt->at(id) < 40.) ){ RecDspT8++;}
				if( (Dspt->at(id) > 40.) and (Dspt->at(id) < 100.) ){ RecDspT9++;}

				if( (abs(Dseta->at(id)) > 0.0) and (abs(Dseta->at(id)) < 0.2) ){ RecDsEta1++;}
				if( (abs(Dseta->at(id)) > 0.2) and (abs(Dseta->at(id)) < 0.4) ){ RecDsEta2++;}
				if( (abs(Dseta->at(id)) > 0.4) and (abs(Dseta->at(id)) < 0.6) ){ RecDsEta3++;}
				if( (abs(Dseta->at(id)) > 0.6) and (abs(Dseta->at(id)) < 0.8) ){ RecDsEta4++;}
				if( (abs(Dseta->at(id)) > 0.8) and (abs(Dseta->at(id)) < 1.0) ){ RecDsEta5++;}
				if( (abs(Dseta->at(id)) > 1.0) and (abs(Dseta->at(id)) < 1.2) ){ RecDsEta6++;}
				if( (abs(Dseta->at(id)) > 1.2) and (abs(Dseta->at(id)) < 1.4) ){ RecDsEta7++;}
				if( (abs(Dseta->at(id)) > 1.4) and (abs(Dseta->at(id)) < 1.6) ){ RecDsEta8++;}
				if( (abs(Dseta->at(id)) > 1.6) and (abs(Dseta->at(id)) < 1.8) ){ RecDsEta9++;}
				if( (abs(Dseta->at(id)) > 1.8) and (abs(Dseta->at(id)) < 2.1) ){ RecDsEta10++;}
			}
			if( deltaR < 0.06){
				hRecEta_Rec2->Fill(abs(Dseta->at(id)));
				hRecPt_Rec2->Fill(Dspt->at(id));						
			}
			if( deltaR < 0.2){
				hRecEta_Rec3->Fill(abs(Dseta->at(id)));
				hRecPt_Rec3->Fill(Dspt->at(id));
			}
			if( deltaR < 0.4){
				hRecEta_Rec4->Fill(abs(Dseta->at(id)));
				hRecPt_Rec4->Fill(Dspt->at(id));			
			}
			t_DsMatching.Fill();		
		}
		}//D* MC protection
		   	
		if (debug)cout << "debug D* 21 --------------------"<< endl;
		//----------------------------------------
		//For D0 MC
		//-----------------------------------------
		if(MCpromptD0mass->size() > 0 && FlagD0MC){
		for(unsigned int i=0; i < MCpromptD0mass->size(); i++)
		{  if (debug)cout << "debug D* 22 --------------------"<< endl;
			/*var_MCpromptD0eta = MCpromptD0eta->at(i);
			var_MCpromptD0phi = MCpromptD0phi->at(i);
			var_MCpromptD0pt = MCpromptD0pt->at(i);
			var_MCpromptD0energy = MCpromptD0energy->at(i);
			var_MCpromptD0p = MCpromptD0p->at(i);
			var_MCpromptD0et = MCpromptD0et->at(i);
			var_MCpromptD0rapidity = MCpromptD0rapidity->at(i);
			var_MCpromptD0mass = MCpromptD0mass->at(i);
			var_MCpromptD0_DispAngle = MCpromptD0_DispAngle->at(i);*/
			var_MCpromptD0lifetime = MCpromptD0lifetime->at(i);

			/*var_MCpromptD0_Keta = MCpromptD0_Keta->at(i);
			var_MCpromptD0_Kphi = MCpromptD0_Kphi->at(i);
			var_MCpromptD0_Kpt = MCpromptD0_Kpt->at(i);
			var_MCpromptD0_Kenergy = MCpromptD0_Kenergy->at(i);
			var_MCpromptD0_Kp = MCpromptD0_Kp->at(i);
			var_MCpromptD0_Ket = MCpromptD0_Ket->at(i);
			var_MCpromptD0_Krapidity = MCpromptD0_Krapidity->at(i);
			var_MCpromptD0_Kmass = MCpromptD0_Kmass->at(i);

			var_MCpromptD0_Pieta = MCpromptD0_Pieta->at(i);
			var_MCpromptD0_Piphi = MCpromptD0_Piphi->at(i);
			var_MCpromptD0_Pipt = MCpromptD0_Pipt->at(i);
			var_MCpromptD0_Pienergy = MCpromptD0_Pienergy->at(i);
			var_MCpromptD0_Pip = MCpromptD0_Pip->at(i);
			var_MCpromptD0_Piet = MCpromptD0_Piet->at(i);
			var_MCpromptD0_Pirapidity = MCpromptD0_Pirapidity->at(i);
			var_MCpromptD0_Pimass = MCpromptD0_Pimass->at(i);*/

			hRecD0Eta_MC->Fill(abs(MCpromptD0eta->at(i)));
			hRecD0Pt_MC->Fill(MCpromptD0pt->at(i));

			if( (MCpromptD0pt->at(i) > 4.) and (MCpromptD0pt->at(i) < 5.) ){ GenD0pT1++;}
			if( (MCpromptD0pt->at(i) > 5.) and (MCpromptD0pt->at(i) < 6.) ){ GenD0pT2++;}
			if( (MCpromptD0pt->at(i) > 6.) and (MCpromptD0pt->at(i) < 7.) ){ GenD0pT3++;}
			if( (MCpromptD0pt->at(i) > 7.) and (MCpromptD0pt->at(i) < 8.) ){ GenD0pT4++;}
			if( (MCpromptD0pt->at(i) > 8.) and (MCpromptD0pt->at(i) < 12.) ){ GenD0pT5++;}
			if( (MCpromptD0pt->at(i) > 12.) and (MCpromptD0pt->at(i) < 16.) ){ GenD0pT6++;}
			if( (MCpromptD0pt->at(i) > 16.) and (MCpromptD0pt->at(i) < 24.) ){ GenD0pT7++;}
			if( (MCpromptD0pt->at(i) > 24.) and (MCpromptD0pt->at(i) < 40.) ){ GenD0pT8++;}
			if( (MCpromptD0pt->at(i) > 40.) and (MCpromptD0pt->at(i) < 100.) ){ GenD0pT9++;}

			if( (abs(MCpromptD0eta->at(i)) > 0.0) and (abs(MCpromptD0eta->at(i)) < 0.2) ){ GenD0Eta1++;}
			if( (abs(MCpromptD0eta->at(i)) > 0.2) and (abs(MCpromptD0eta->at(i)) < 0.4) ){ GenD0Eta2++;}
			if( (abs(MCpromptD0eta->at(i)) > 0.4) and (abs(MCpromptD0eta->at(i)) < 0.6) ){ GenD0Eta3++;}
			if( (abs(MCpromptD0eta->at(i)) > 0.6) and (abs(MCpromptD0eta->at(i)) < 0.8) ){ GenD0Eta4++;}
			if( (abs(MCpromptD0eta->at(i)) > 0.8) and (abs(MCpromptD0eta->at(i)) < 1.0) ){ GenD0Eta5++;}
			if( (abs(MCpromptD0eta->at(i)) > 1.0) and (abs(MCpromptD0eta->at(i)) < 1.2) ){ GenD0Eta6++;}
			if( (abs(MCpromptD0eta->at(i)) > 1.2) and (abs(MCpromptD0eta->at(i)) < 1.4) ){ GenD0Eta7++;}
			if( (abs(MCpromptD0eta->at(i)) > 1.4) and (abs(MCpromptD0eta->at(i)) < 1.6) ){ GenD0Eta8++;}
			if( (abs(MCpromptD0eta->at(i)) > 1.6) and (abs(MCpromptD0eta->at(i)) < 1.8) ){ GenD0Eta9++;}
			if( (abs(MCpromptD0eta->at(i)) > 1.8) and (abs(MCpromptD0eta->at(i)) < 2.1) ){ GenD0Eta10++;}

			//t_D0MC.Fill();
			if (debug)cout << "debug D* 23 --------------------"<< endl;
			//----------------------------------------
			//For D0 Matching For Recosntruction Eff
			//----------------------------------------
			if ( D0Kpimass->size() == 0) continue;
			double minDeltaR = 1000.;
			int id = -1;
			double deltaR = 0.;
			for(unsigned int j = 0; j < D0Kpimass->size(); j++)
			{  													
				if( D0Kpi_VtxProb->at(j) < 0.01) continue;
				if( D0KpiDispAngle->at(j) < 0.99 ) continue;
				if( abs(D0Kpimass->at(j)-1.86484) > 0.1) continue;
				if( D0Kpipt->at(j) < 2. ) continue;	
				if( D0Kpis3D->at(j) < SigCutD0) continue;
				
				//DeltaR Function
				deltaR = sqrt( pow( (MCpromptD0eta->at(i)-D0Kpieta->at(j)) ,2) + pow( (MCpromptD0phi->at(i)-D0Kpiphi->at(j)) ,2) );
				if( deltaR < minDeltaR) { id = j; minDeltaR = deltaR;} //Get the lower delta
				if (debug)cout << "debug D0 Matching --------------------"<< endl;	
			}	

			if (id == -1) continue;
			
			var_MCpromptD0lifetimeMatching = MCpromptD0lifetime->at(i);
			/*var_D0etaMatching = D0etaMatching->at(id);
			var_D0phiMatching = D0phiMatching->at(id);
			var_D0ptMatching = D0ptMatching->at(id);
			var_D0KtMatching = D0KtMatching->at(id);
			var_D0SxyMatching = D0SxyMatching->at(id);
			var_D0OpAngleMatching = D0OpAngleMatching->at(id);

			var_MCD0etaMatching = MCD0etaMatching->at(id);
			var_MCD0phiMatching = MCD0phiMatching->at(id);
			var_MCD0ptMatching = MCD0ptMatching->at(id);*/
			var_deltaRD0Matching = deltaR;
		
			D0deltaEta = MCpromptD0eta->at(i)-D0Kpieta->at(id);
			D0deltaPhi = MCpromptD0phi->at(i)-D0Kpiphi->at(id);
			D0deltaPt = MCpromptD0pt->at(i)-D0Kpipt->at(id);

			double recMass = -1; double GenMass = -1; // Protection
			double recPt = -1; double GenPt = -1;
			double recEta = -1; double GenEta = -1;
			double recTime = -1; double GenTime = -1;
								
			if( deltaR < 0.03)
			{	hRecD0Eta_Rec1->Fill( abs(D0Kpieta->at(id)) );
				hRecD0Pt_Rec1->Fill( D0Kpipt->at(id) );

				recMass = D0Kpimass->at(id); GenMass = MCpromptD0mass->at(i);
				recPt = D0Kpipt->at(id); GenPt = MCpromptD0pt->at(i);
				recEta = D0Kpieta->at(id); GenEta = MCpromptD0eta->at(i);
				recTime = (D0Kpimass->at(id)*D0KpidXY->at(id)*pow(10,-2)) / (c*D0Kpipt->at(id)); 
				GenTime = MCpromptD0lifetime->at(i);

				hdeltaRD0->Fill( recMass, GenMass );
				hdeltaRD0Pt->Fill( recPt, GenPt );
				hdeltaRD0Eta->Fill( abs(recEta), abs(GenEta) );
				hdeltaRD0Time->Fill( recTime, GenTime );

				if( (D0Kpipt->at(id) > 4.) and (D0Kpipt->at(id) < 5.) ){ RecD0pT1++;}
				if( (D0Kpipt->at(id) > 5.) and (D0Kpipt->at(id) < 6.) ){ RecD0pT2++;}
				if( (D0Kpipt->at(id) > 6.) and (D0Kpipt->at(id) < 7.) ){ RecD0pT3++;}
				if( (D0Kpipt->at(id) > 7.) and (D0Kpipt->at(id) < 8.) ){ RecD0pT4++;}
				if( (D0Kpipt->at(id) > 8.) and (D0Kpipt->at(id) < 12.) ){ RecD0pT5++;}
				if( (D0Kpipt->at(id) > 12.) and (D0Kpipt->at(id) < 16.) ){ RecD0pT6++;}
				if( (D0Kpipt->at(id) > 16.) and (D0Kpipt->at(id) < 24.) ){ RecD0pT7++;}
				if( (D0Kpipt->at(id) > 24.) and (D0Kpipt->at(id) < 40.) ){ RecD0pT8++;}
				if( (D0Kpipt->at(id) > 40.) and (D0Kpipt->at(id) < 100.) ){ RecD0pT9++;}

				if( (abs(D0Kpieta->at(id)) > 0.0) and (abs(D0Kpieta->at(id)) < 0.2) ){ RecD0Eta1++;}
				if( (abs(D0Kpieta->at(id)) > 0.2) and (abs(D0Kpieta->at(id)) < 0.4) ){ RecD0Eta2++;}
				if( (abs(D0Kpieta->at(id)) > 0.4) and (abs(D0Kpieta->at(id)) < 0.6) ){ RecD0Eta3++;}
				if( (abs(D0Kpieta->at(id)) > 0.6) and (abs(D0Kpieta->at(id)) < 0.8) ){ RecD0Eta4++;}
				if( (abs(D0Kpieta->at(id)) > 0.8) and (abs(D0Kpieta->at(id)) < 1.0) ){ RecD0Eta5++;}
				if( (abs(D0Kpieta->at(id)) > 1.0) and (abs(D0Kpieta->at(id)) < 1.2) ){ RecD0Eta6++;}
				if( (abs(D0Kpieta->at(id)) > 1.2) and (abs(D0Kpieta->at(id)) < 1.4) ){ RecD0Eta7++;}
				if( (abs(D0Kpieta->at(id)) > 1.4) and (abs(D0Kpieta->at(id)) < 1.6) ){ RecD0Eta8++;}
				if( (abs(D0Kpieta->at(id)) > 1.6) and (abs(D0Kpieta->at(id)) < 1.8) ){ RecD0Eta9++;}
				if( (abs(D0Kpieta->at(id)) > 1.8) and (abs(D0Kpieta->at(id)) < 2.1) ){ RecD0Eta10++;}

			}
			if( deltaR < 0.06)
			{	hRecD0Eta_Rec2->Fill( abs(D0Kpieta->at(id)) );
				hRecD0Pt_Rec2->Fill( D0Kpipt->at(id) );
			
			}
			if( deltaR < 0.2)
			{	hRecD0Eta_Rec3->Fill( abs(D0Kpieta->at(id)) );
				hRecD0Pt_Rec3->Fill( D0Kpipt->at(id) );
			}

			if( deltaR < 0.4)
			{	hRecD0Eta_Rec4->Fill( abs(D0Kpieta->at(id)) );
				hRecD0Pt_Rec4->Fill( D0Kpipt->at(id) );
			}
		
			t_D0Matching.Fill(); //Fill tree
		}
		}			
#endif
		if (debug)cout << "debug 24 --------------------"<< endl;	
	}//End loop tree entries for file f1

	if (debug)cout << "debug 25 --------------------"<< endl;

	cout << "--------------------"<< endl;
	cout << "NumberOfEvents: " << NumberOfEvents << endl;
	cout << "EventsAfterTrigger: " << EventsAfterTrigger << endl;
	cout << "--------------------"<< endl;
	cout << "PossibleDs: " << PossibleDs << endl;
	cout << "DsD0ProbVtx: " << DsD0ProbVtx << endl;
	cout << "DsD0Angle: " << DsD0Angle << endl;
	cout << "D0fromDsMinusPDG: " << D0fromDsMinusPDG << endl;
	cout << "DsD0pt: " << DsD0pt << endl;
	cout << "DsD0deltaM: " << DsD0deltaM << endl;
	cout << "DsD0Sig: " << DsD0Sig << endl;
	cout << "--------------------"<< endl;
	cout << "PossibleD0: " << PossibleD0 <<  endl;
	cout << "D0ProbVtx: " << D0ProbVtx <<  endl;
	cout << "D0Angle: " << D0Angle << endl;
	cout << "CountD0minusPDG: " << CountD0minusPDG << endl;
	cout << "CountD0pt: " << CountD0pt  << endl;
	cout << "D0Sig: " << D0Sig << endl;
	
	if (Dataset == "MC_DStarToD0Pi_D0KPi" or Dataset == "MC_MinBias"){
	cout << "------Rec Ds Efficiency--------------"<< endl;
	cout << "pT range / N Gen / N Rec / Efficiency  ----------"<< endl;
	cout << "$4 < pT < 5$ & " << GenDspT1 << " & "<< RecDspT1 << " & " << RecDspT1*1./GenDspT1 << " end" << endl;
	cout << "$5 < pT < 6$ & " << GenDspT2 << " & "<< RecDspT2 << " & " << RecDspT2*1./GenDspT2 << " end" << endl;
	cout << "$6 < pT < 7$ & " << GenDspT3 << " & "<< RecDspT3 << " & " << RecDspT3*1./GenDspT3 << " end" << endl;
	cout << "$7 < pT < 8$ & " << GenDspT4 << " & "<< RecDspT4 << " & " << RecDspT4*1./GenDspT4 << " end" << endl;
	cout << "$8 < pT < 12$ & " << GenDspT5 << " & "<< RecDspT5 << " & " << RecDspT5*1./GenDspT5 << " end" << endl;
	cout << "$12 < pT < 16$ & " << GenDspT6 << " & "<< RecDspT6 << " & " << RecDspT6*1./GenDspT6 << " end" << endl;
	cout << "$16 < pT < 24$ & " << GenDspT7 << " & "<< RecDspT7 << " & " << RecDspT7*1./GenDspT7 << " end" << endl;
	cout << "$24 < pT < 40$ & " << GenDspT8 << " & "<< RecDspT8 << " & " << RecDspT8*1./GenDspT8 << " end" << endl;
	cout << "$40 < pT < 100$ & " << GenDspT9 << " & "<< RecDspT9 << " & " << RecDspT9*1./GenDspT9 << " end" << endl;
	cout << "eta range / N Gen / N Rec / Efficiency  ----------"<< endl;
	cout << "$0.0 < eta < 0.2 $ & " << GenDsEta1 << " & "<< RecDsEta1 << " & " << RecDsEta1*1./GenDsEta1 << " end" << endl;
	cout << "$0.2 < eta < 0.4 $ & " << GenDsEta2 << " & "<< RecDsEta2 << " & " << RecDsEta2*1./GenDsEta2 << " end" << endl;
	cout << "$0.4 < eta < 0.6 $ & " << GenDsEta3 << " & "<< RecDsEta3 << " & " << RecDsEta3*1./GenDsEta3 << " end" << endl;
	cout << "$0.6 < eta < 0.7 $ & " << GenDsEta4 << " & "<< RecDsEta4 << " & " << RecDsEta4*1./GenDsEta4 << " end" << endl;
	cout << "$0.8 < eta < 1.0 $ & " << GenDsEta5 << " & "<< RecDsEta5 << " & " << RecDsEta5*1./GenDsEta5 << " end" << endl;
	cout << "$1.0 < eta < 1.2 $ & " << GenDsEta6 << " & "<< RecDsEta6 << " & " << RecDsEta6*1./GenDsEta6 << " end" << endl;
	cout << "$1.2 < eta < 1.4 $ & " << GenDsEta7 << " & "<< RecDsEta7 << " & " << RecDsEta7*1./GenDsEta7 << " end" << endl;
	cout << "$1.4 < eta < 1.6 $ & " << GenDsEta8 << " & "<< RecDsEta8 << " & " << RecDsEta8*1./GenDsEta8 << " end" << endl;
	cout << "$1.6 < eta < 1.8 $ & " << GenDsEta9 << " & "<< RecDsEta9 << " & " << RecDsEta9*1./GenDsEta9 << " end" << endl;
	cout << "$1.8 < eta < 2.1 $ & " << GenDsEta10 << " & "<< RecDsEta10 << " & " << RecDsEta10*1./GenDsEta10 << " end" << endl;
	cout << "------Rec D0 Efficiency--------------"<< endl;
	cout << "pt range / N Gen / N Rec / Efficiency  ----------"<< endl;
	cout << "pT range / N Gen / N Rec / Efficiency  ----------"<< endl;
	cout << "$4 < pT < 5$ & " << GenD0pT1 << " & "<< RecD0pT1 << " & " << RecD0pT1*1./GenD0pT1 << " end" << endl;
	cout << "$5 < pT < 6$ & " << GenD0pT2 << " & "<< RecD0pT2 << " & " << RecD0pT2*1./GenD0pT2 << " end" << endl;
	cout << "$6 < pT < 7$ & " << GenD0pT3 << " & "<< RecD0pT3 << " & " << RecD0pT3*1./GenD0pT3 << " end" << endl;
	cout << "$7 < pT < 8$ & " << GenD0pT4 << " & "<< RecD0pT4 << " & " << RecD0pT4*1./GenD0pT4 << " end" << endl;
	cout << "$8 < pT < 12$ & " << GenD0pT5 << " & "<< RecD0pT5 << " & " << RecD0pT5*1./GenD0pT5 << " end" << endl;
	cout << "$12 < pT < 16$ & " << GenD0pT6 << " & "<< RecD0pT6 << " & " << RecD0pT6*1./GenD0pT6 << " end" << endl;
	cout << "$16 < pT < 24$ & " << GenD0pT7 << " & "<< RecD0pT7 << " & " << RecD0pT7*1./GenD0pT7 << " end" << endl;
	cout << "$24 < pT < 40$ & " << GenD0pT8 << " & "<< RecD0pT8 << " & " << RecD0pT8*1./GenD0pT8 << " end" << endl;
	cout << "$40 < pT < 100$ & " << GenD0pT9 << " & "<< RecD0pT9 << " & " << RecD0pT9*1./GenD0pT9 << " end" << endl;
	cout << "eta range / N Gen / N Rec / Efficiency  ----------"<< endl;
	cout << "$0.0 < eta < 0.2 $ & " << GenD0Eta1 << " & "<< RecD0Eta1 << " & " << RecD0Eta1*1./GenD0Eta1 << " end" << endl;
	cout << "$0.2 < eta < 0.4 $ & " << GenD0Eta2 << " & "<< RecD0Eta2 << " & " << RecD0Eta2*1./GenD0Eta2 << " end" << endl;
	cout << "$0.4 < eta < 0.6 $ & " << GenD0Eta3 << " & "<< RecD0Eta3 << " & " << RecD0Eta3*1./GenD0Eta3 << " end" << endl;
	cout << "$0.6 < eta < 0.7 $ & " << GenD0Eta4 << " & "<< RecD0Eta4 << " & " << RecD0Eta4*1./GenD0Eta4 << " end" << endl;
	cout << "$0.8 < eta < 1.0 $ & " << GenD0Eta5 << " & "<< RecD0Eta5 << " & " << RecD0Eta5*1./GenD0Eta5 << " end" << endl;
	cout << "$1.0 < eta < 1.2 $ & " << GenD0Eta6 << " & "<< RecD0Eta6 << " & " << RecD0Eta6*1./GenD0Eta6 << " end" << endl;
	cout << "$1.2 < eta < 1.4 $ & " << GenD0Eta7 << " & "<< RecD0Eta7 << " & " << RecD0Eta7*1./GenD0Eta7 << " end" << endl;
	cout << "$1.4 < eta < 1.6 $ & " << GenD0Eta8 << " & "<< RecD0Eta8 << " & " << RecD0Eta8*1./GenD0Eta8 << " end" << endl;
	cout << "$1.6 < eta < 1.8 $ & " << GenD0Eta9 << " & "<< RecD0Eta9 << " & " << RecD0Eta9*1./GenD0Eta9 << " end" << endl;
	cout << "$1.8 < eta < 2.1 $ & " << GenD0Eta10 << " & "<< RecD0Eta10 << " & " << RecD0Eta10*1./GenD0Eta10 << " end" << endl;
	cout << "--------------------"<< endl;

	cout << "----------  C O N T A M I N A T I O N   D * -------------"<< endl;
	cout << "Total / Dstar Prompt / Dstar Secondary / Contamination  ----------"<< endl;
	cout << DstarTotal << " & "<< DstarTotal-DstarSec << " & "<< DstarSec  << " & " << DstarSec/DstarTotal << endl;
	}
	

	TCanvas* canvas = new TCanvas("PreliminarStudyDs","",900,600);
	canvas->Divide(3,3);
	canvas->cd(1);
	D0VtxProbHisto->Draw("Histo");
	canvas->cd(2);
	TotalD0fromDSs3DHisto->Draw("Histo");
	canvas->cd(3);
	D0fromDSs3DHisto->Draw("Histo");	
	canvas->cd(4);
	DsmassHisto->Draw("e1p");
	canvas->cd(5);
	D0massHisto->Draw("e1p");
	canvas->cd(6);
	DsMinusD0Histo->Draw("Histo");
	canvas->cd(7);
	DslifetimeHisto->Draw("Histo");
	canvas->cd(8);
	canvas->SaveAs(Form("%s%s/PreliminarStudyDs_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));

	canvas = new TCanvas("PreliminarStudyD0","",900,600);
	canvas->Divide(3,3);
	canvas->cd(1);
	D0Kpi_VtxProbHisto->Draw("Histo");
	canvas->cd(2);
	TotalD0Kpis3DHistoHisto->Draw("Histo");
	canvas->cd(3);
	D0Kpis3DHisto->Draw("Histo");
	canvas->cd(4);	
	D0KpimassHisto->Draw("e1p");
	canvas->cd(5);
	D0lifetimeHisto->Draw("Histo");
	canvas->cd(6);
	canvas->cd(7);
	canvas->cd(8);
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));

	canvas = new TCanvas("PreliminarStudyDsMass","",900,600);
	canvas->Divide(3,3);
	canvas->cd(1);
	DsmassHisto0->Draw("e1p");
	TLatex* tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 0.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig = 0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig = 0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(2);
	DsmassHisto1->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 0.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 0.5) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 0.5) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(3);
	DsmassHisto2->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 1.0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 1.0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 1.0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(4);	
	DsmassHisto3->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 1.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 1.5) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 1.5) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(5);
	DsmassHisto4->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 2.0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 2.0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 2.0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(6);
	DsmassHisto5->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 2.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 2.5) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 2.5) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(7);
	DsmassHisto6->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 3.0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 3.0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 3.0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(8);
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMass_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));

	canvas = new TCanvas("PreliminarStudyDsMassPt","",900,600);
	canvas->Divide(3,3);
	canvas->cd(1);
	DsmassPtHisto1->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT = 0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT = 0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(2);
	DsmassPtHisto2->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 0.5 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 0.5 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(3);
	DsmassPtHisto3->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 1.0 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 1.0 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(4);
	DsmassPtHisto4->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 1.5 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 1.5 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(5);
	DsmassPtHisto5->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 2.0 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 2.0 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(6);
	DsmassPtHisto6->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 2.5 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 2.5 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(7);
	DsmassPtHisto7->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 3.0 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 3.0 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(8);
	DsmassPtHisto8->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 3.5 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 3.5 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(9);
	DsmassPtHisto9->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 4.0 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 4.0 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMassPt_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMassPt_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));

	canvas = new TCanvas("PreliminarStudyD0Mass","",900,600);
	canvas->Divide(3,4);
	canvas->cd(1);
	D0KpimassHisto0->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig = 0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig = 0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig = 0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(2);
	D0KpimassHisto1->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 0.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 0.5) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 0.5) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(3);
	D0KpimassHisto2->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 1.0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 1.0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 1.0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(4);	
	D0KpimassHisto3->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 1.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 1.5) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 1.5) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(5);
	D0KpimassHisto4->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 2.0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 2.0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 2.0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(6);
	D0KpimassHisto5->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 2.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 2.5) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 2.5) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(7);
	D0KpimassHisto6->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 3.0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 3.0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 3.0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(8);
	D0KpimassHisto7->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 3.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 3.5) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 3.5) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(9);
	D0KpimassHisto8->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 4.0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 4.0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 4.0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(10);
	D0KpimassHisto9->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 4.5) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 4.5) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 4.5) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(11);
	D0KpimassHisto10->Draw("e1p");
	tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig > 5.0) }");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig > 5.0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig > 5.0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0Mass_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0Mass_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));

	if (debug)cout << "debug 15 --------------------"<< endl;

	canvas = new TCanvas("PreliminarStudyD0MassPt","",900,600);
	canvas->Divide(3,3);
	canvas->cd(1);
	D0KpimassPtHisto1->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT = 0) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT = 0) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(2);
	D0KpimassPtHisto2->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 0.5 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 0.5 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(3);
	D0KpimassPtHisto3->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 1.0 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 1.0 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(4);
	D0KpimassPtHisto4->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 1.5 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 1.5 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(5);
	D0KpimassPtHisto5->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 2.0 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 2.0 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(6);
	D0KpimassPtHisto6->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 2.5 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 2.5 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(7);
	D0KpimassPtHisto7->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 3.0 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 3.0 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(8);
	D0KpimassPtHisto8->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 3.5 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 3.5 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->cd(9);
	D0KpimassPtHisto9->Draw("e1p");
	if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT > 4.0 GeV) }");}
	if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT > 4.0 GeV) }");}
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0MassPt_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0MassPt_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));

//----------------------------------------------------------------------------------------------------------------
#if GENvaraibles == 0
#elif GENvaraibles == 1
	canvas = new TCanvas("StudyRecGenDstar","",900,600);
	canvas->Divide(2,2);
	//--------
	canvas->cd(1);	hdeltaRDstar->Draw();
	tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	//--------
	canvas->cd(2);	hdeltaRDstarPt->Draw();
	tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	//--------
	canvas->cd(3);	hdeltaRDstarEta->Draw();
	tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	//--------
	canvas->cd(4);	hdeltaRDstarTime->Draw();
	tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	//--------
	canvas->SaveAs(Form("%s%s/StudyRecGenDstar_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/StudyRecGenDstar_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	//**********************************************************
	canvas = new TCanvas("StudyRecGenD0","",900,600);
	canvas->Divide(2,2);
	//--------
	canvas->cd(1);	hdeltaRD0->Draw();
	tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	//--------
	canvas->cd(2);	hdeltaRD0Pt->Draw();
	tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	//--------
	canvas->cd(3);	hdeltaRD0Eta->Draw();
	tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	//--------
	canvas->cd(4);	hdeltaRD0Time->Draw();
	tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	//--------
	canvas->SaveAs(Form("%s%s/StudyRecGenD0_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/StudyRecGenD0_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
#endif

	//t_analysis.Branch("D0mass",&D0mass);
	//TBranch *D0mass_branch; D0mass_branch = 
	f_analysis.cd();

	t_analysis.Write();  //Write in the root file
	t_DstarWrongCombination.Write();
	t_D0analysis.Write();

	DsmassHisto1->Write();
	DsmassHisto2->Write();
	DsmassHisto3->Write();
	DsmassHisto4->Write();
	DsmassHisto5->Write(); 
	DsmassHisto6->Write();

	DsmassPtHisto1->Write();
	DsmassPtHisto2->Write();
	DsmassPtHisto3->Write();
	DsmassPtHisto4->Write();
	DsmassPtHisto5->Write(); 
	DsmassPtHisto6->Write();
	DsmassPtHisto7->Write();
	DsmassPtHisto8->Write();
	DsmassPtHisto9->Write();
		
	D0KpimassHisto1->Write();
	D0KpimassHisto2->Write();
	D0KpimassHisto3->Write();
	D0KpimassHisto4->Write();
	D0KpimassHisto5->Write();
	D0KpimassHisto6->Write();
	D0KpimassHisto7->Write();
	D0KpimassHisto8->Write();
	D0KpimassHisto9->Write();
	D0KpimassHisto10->Write();

	D0KpimassPtHisto1->Write();
	D0KpimassPtHisto2->Write();
	D0KpimassPtHisto3->Write();
	D0KpimassPtHisto4->Write();
	D0KpimassPtHisto5->Write();
	D0KpimassPtHisto6->Write();
	D0KpimassPtHisto7->Write();
	D0KpimassPtHisto8->Write();
	D0KpimassPtHisto9->Write();
	
#if GENvaraibles == 0	
#elif GENvaraibles == 1
	t_DsMC.Write(); //Write() if a file is open, this function writes a root objectics on it.
	t_D0MC.Write();
	t_DsMatching.Write();
	t_D0Matching.Write();

	hdeltaRDstar->Write();
	hdeltaRDstarPt->Write();
	hdeltaRDstarEta->Write();
	hdeltaRDstarTime->Write();
	hdeltaRD0->Write();
	hdeltaRD0Pt->Write();
	hdeltaRD0Eta->Write();
	hdeltaRD0Time->Write();

	hRecEta_MC->Write(); 
	hRecPt_MC->Write();
	hRecEta_Rec1->Write();
	hRecPt_Rec1->Write();
	hRecEta_Rec2->Write();
	hRecPt_Rec2->Write();
	hRecEta_Rec3->Write();
	hRecPt_Rec3->Write();
	hRecEta_Rec4->Write();
	hRecPt_Rec4->Write();
	//------------------------------------------
	hRecD0Eta_MC->Write();
	hRecD0Pt_MC->Write();
	hRecD0Eta_Rec1->Write();
	hRecD0Pt_Rec1->Write();
	hRecD0Eta_Rec2->Write();
	hRecD0Pt_Rec2->Write();
	hRecD0Eta_Rec3->Write();
	hRecD0Pt_Rec3->Write();
	hRecD0Eta_Rec4->Write();
	hRecD0Pt_Rec4->Write();
#endif	
	//------------------------------------------
	f_Puanalysis.cd();

	//D* Pileup weight
	hDsmass->Write();
	hDseta->Write();
	hDsphi->Write();
	hDspt->Write();
	hDslifetime->Write();

	hDsMinusD0->Write();
	hD0mass->Write();
	hD0eta->Write();
	hD0phi->Write();
	hD0pt->Write();

	hTrkKpt->Write();
	hTrkKchi2->Write();
	hTrkKnhits->Write();
	hTrkKdz->Write();
	hTrkKdxy->Write();
	hTrkKeta->Write();
	hTrkKphi->Write();

	hTrkSpt->Write();
	hTrkSchi2->Write();
	hTrkSnhits->Write();
	hTrkSdz->Write();
	hTrkSdxy->Write();
	hTrkSeta->Write();
	hTrkSphi->Write();

	hTrkpipt->Write();
	hTrkpichi2->Write();
	hTrkpinhits->Write();
	hTrkpidz->Write();
	hTrkpidxy->Write();
	hTrkpieta->Write();
	hTrkpiphi->Write();

	//D0 Pileup weight
	hD0Kpimass->Write();
	hD0Kpieta->Write();
	hD0Kpiphi->Write();
	hD0Kpipt->Write();
	hD0lifetime->Write();

	hTrkD0Kpt->Write();
	hTrkD0Kchi2->Write();
	hTrkD0Knhits->Write();
	hTrkD0Kdz->Write();
	hTrkD0Kdxy->Write();
	hTrkD0Keta->Write();
	hTrkD0kphi->Write();

	hTrkD0pipt->Write();
	hTrkD0pichi2->Write();
	hTrkD0pinhits->Write();
	hTrkD0pidz->Write();
	hTrkD0pidxy->Write();
	hTrkD0pieta->Write();
	hTrkD0piphi->Write();

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	printf("Time taken: %.2fm\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC*60.));
	printf("Time taken: %.2fh\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC*60.*60));
	cout << "-------END PROGRAM-------------"<< endl;

}//end program
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
	const char *ext =".root"; //extenson of files you want add	
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
			{	//cout << "File: "<< fname << endl; 
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
	hh->SetTitle(TITLES); hh->SetMarkerStyle(7);
	hh->SetFillColor(COLOR); hh->Sumw2();
	//hh->Draw(TYPE); const char* TYPE
   return hh ;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH1* makeTH1Rec(const char* name, const int nbins_low , double* nbins_high ){
	
	//Create a TH1D (one dimension) with alternatives bins
	TH1D *hh = new TH1D(name, name, nbins_low, nbins_high);
	hh->SetMarkerStyle(7); hh->Sumw2(); return hh ;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//TBranch* makeBranch(TTree *myTree, const char* name, vector<double>* VECTOR) 
//{
//	TBranch *branch = myTree->GetBranch(name);
//	branch->SetAddress(&VECTOR);
//	return branch ;
//}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
