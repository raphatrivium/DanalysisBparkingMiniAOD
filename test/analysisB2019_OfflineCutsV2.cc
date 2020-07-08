/////////////////////////////////////////////////////////////////
//
// Program to extract data (variables and vector format)
// of the root file and create histograms. This root file
// is made of others root files, therefore, has multiple
// entries which we access with a diferent method than a
// normal root file.
//
// by Raphael Gomes de Souza - UERJ - RJ - Brazil
////////////////////////////////////////////////////////////
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
#include "TText.h"
#include "TAttLine.h"
//#ifndef __CINT__
#include "TCanvas.h"
#include "TAxis.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TChain.h"
#include <time.h>
#include "TLegend.h"
#include <TLatex.h>
#include <typeinfo>
#include "TEfficiency.h"

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
TH2* makeTH2GenRec( TString name, TString nameAxis, Double_t BIN, Double_t NMIN , Double_t NMAX );
void TGraphAsymmErrorsRec( TH1 *h1, TH1 *h2, const char* hTITLE, string path);
//--------------------------------------------------------------
// M A I N   F U N C T I O N 
//--------------------------------------------------------------
void analysisB2019_OfflineCutsV2()
{
	clock_t tStart = clock();

	//Datasets: BparkingDataset1 - MC_DStarToD0Pi_D0KPi - MC_MinBias
	string Dataset = "MC_DStarToD0Pi_D0KPi";
	//Put "" for select no trigger 
	string TriggerPath = ""; //HLT_Mu9_IP6
	if (Dataset != "MC_DStarToD0Pi_D0KPi" and Dataset != "MC_MinBias"){TriggerPath = "HLT_Mu9_IP6";} //HLT_Mu9_IP6
	//Save the file
	string path = "/eos/user/r/ragomesd/crab/haddteste/";
	string path2 = "/eos/user/r/ragomesd/analysisB2019/canvas/";
	double SigCutD0fromDs, SigCutD0;

	if (Dataset != "MC_DStarToD0Pi_D0KPi"){SigCutD0fromDs = 3.; SigCutD0 = 5.; } //3
	if (Dataset == "MC_DStarToD0Pi_D0KPi" or Dataset == "MC_MinBias"){SigCutD0fromDs = 3.; SigCutD0 = 5.; }

	bool FlagDs = true;
	bool FlagDsMC = true;
	bool FlagD0 = true;
	bool FlagD0MC = true;
	bool debug = false;
		
	//TChain* t1 = callTchain("/eos/user/r/ragomesd/crab/ParkingBPH2/ParkingBPH_Run2018A_MINIAOD/200105_042547/0000/");
	TChain* t1 = callTchain("./");

	//**********************************************************		
	Long64_t nentries = t1->GetEntries(); //Reading Number of tree entries of the file
	nentries = 200000; //Test
	cout<< "Number of tree entries: "<< nentries << endl;
	Long64_t partpercent = nentries*0.05; //Percent done to show
	double c = 299792458; //Light speed [m/s]
	//-------Create root file-------------------------------------
	//TFile f_analysis("DMesonsHistograms.root","recreate"); //Creates root file
	TFile f_analysis(Form("%s%s.root", path.c_str(), Dataset.c_str() ),"recreate"); //Creates root file
	TTree t_analysis("t_analysis","t_analysis"); //Creates a Tree
	TTree t_DstarWrongCombination("t_DstarWrongCombination","t_DstarWrongCombination");
	TTree t_D0analysis("t_D0analysis","t_D0analysis");
#if GENvaraibles == 0
#elif GENvaraibles == 1
	TTree t_DsMC("t_DsMC","t_DsMC"); 
	TTree t_D0MC("t_D0MC","t_D0MC");
	TTree t_DsMatching("t_DsMatching","t_DsMatching");
	TTree t_D0Matching("t_D0Matching","t_D0Matching");
	TTree t_DstarContamination("t_DstarContamination","t_DstarContamination");
	TTree t_D0Contamination("t_D0Contamination","t_D0Contamination");
#endif
	int sizeArray = 0;
	
	//----------------------------------------
	//For D*
	//-----------------------------------------
	TString DstarString [] = {"Dsmass", "Dseta", "Dsphi", "Dspt", "DSDeltaR", "D0mass", "D0eta", "D0phi", "D0pt", "D0fromDSdXY", "D0fromDSsXY", "D0fromDSs3D", "D0fromDSd3D", "Anglephi", "D0_VtxProb"};
	//----------------------------------------
	// Vectors must be equal to "0", not "0.". Possible segmental error if you compile it with c++.
	vector< vector<double>* > DstarVec;
	vector<double>* Dsmass =0; vector<double>* Dseta =0; vector<double>* Dsphi =0; vector<double>* Dspt =0; vector<double>* DSDeltaR =0; vector<double>* D0mass =0; vector<double>* D0eta =0; vector<double>* D0phi =0; vector<double>* D0pt =0; vector<double>* D0fromDSdXY =0; vector<double>* D0fromDSsXY =0; vector<double>* D0fromDSs3D =0; vector<double>* D0fromDSd3D =0; vector<double>* Anglephi =0; vector<double>* D0_VtxProb =0;
	DstarVec.push_back(Dsmass); DstarVec.push_back(Dseta); DstarVec.push_back(Dsphi); DstarVec.push_back(Dspt); DstarVec.push_back(DSDeltaR); DstarVec.push_back(D0mass); DstarVec.push_back(D0eta); DstarVec.push_back(D0phi); DstarVec.push_back(D0pt); DstarVec.push_back(D0fromDSdXY); DstarVec.push_back(D0fromDSsXY); DstarVec.push_back(D0fromDSs3D); DstarVec.push_back(D0fromDSd3D); DstarVec.push_back(Anglephi); DstarVec.push_back(D0_VtxProb);
	// Input Branches
	vector< TBranch* > DstarBranch; DstarBranch.clear();
	TBranch *b_Dsmass; TBranch *b_Dseta; TBranch *b_Dsphi; TBranch *b_Dspt; TBranch *b_DSDeltaR; TBranch *b_D0mass; TBranch *b_D0eta; TBranch *b_D0phi; TBranch *b_D0pt; TBranch *b_D0fromDSdXY; TBranch *b_D0fromDSsXY; TBranch *b_D0fromDSs3D; TBranch *b_D0fromDSd3D; TBranch *b_Anglephi; TBranch *b_D0_VtxProb;
	DstarBranch.push_back(b_Dsmass); DstarBranch.push_back(b_Dseta); DstarBranch.push_back(b_Dsphi); DstarBranch.push_back(b_Dspt); DstarBranch.push_back(b_DSDeltaR); DstarBranch.push_back(b_D0mass); DstarBranch.push_back(b_D0eta); DstarBranch.push_back(b_D0phi); DstarBranch.push_back(b_D0pt); DstarBranch.push_back(b_D0fromDSdXY); DstarBranch.push_back(b_D0fromDSsXY); DstarBranch.push_back(b_D0fromDSs3D); DstarBranch.push_back(b_D0fromDSd3D); DstarBranch.push_back(b_Anglephi); DstarBranch.push_back(b_D0_VtxProb); 		
	// Works but has a warning message: "warning: ‘object’ may be used uninitialized in this function"
	//TBranch* DstarBranchvec[] = {b_Dsmass, b_D0mass};
	//int sizeArray = sizeof(DstarBranchvec)/sizeof(DstarBranchvec[0]); // Number of elements
	//vector< TBranch* > DstarBranch;  DstarBranch.clear();
	//DstarBranch.insert(DstarBranch.begin(), DstarBranchvec, DstarBranchvec+sizeArray );
	//----------------------------------------------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < DstarBranch.size(); k++)
	{ t1->SetBranchAddress( DstarString[k], &DstarVec[k], &DstarBranch[k]); }
	//----------------------------------------------------------------------------
	//Variables 
	double var_Dsmass =0.; double var_Dseta =0.; double var_Dsphi =0.; double var_Dspt =0.; double var_DSDeltaR =0.; double var_D0mass =0.; double var_D0eta =0.; double var_D0phi =0.; double var_D0pt =0.; double var_D0fromDSdXY =0.; double var_D0fromDSsXY =0.; double var_D0fromDSs3D =0.; double var_D0fromDSd3D =0.; double var_Anglephi =0.; double var_D0_VtxProb =0.;  
	double DstarVar[] = {var_Dsmass, var_Dseta, var_Dsphi, var_Dspt, var_DSDeltaR, var_D0mass, var_D0eta, var_D0phi, var_D0pt, var_D0fromDSdXY, var_D0fromDSsXY, var_D0fromDSs3D, var_D0fromDSd3D, var_Anglephi, var_D0_VtxProb, };
	sizeArray = sizeof(DstarVar)/sizeof(DstarVar[0]); // Number of elements
	vector<double> DstarVarVec; DstarVarVec.clear();
	DstarVarVec.insert(DstarVarVec.begin(), DstarVar, DstarVar+sizeArray );
	//----------------------------------------------------------------------------
	// CREATING OUTPUT BRANCHES
	for( unsigned int k = 0; k < DstarVarVec.size(); k++)
	{ t_analysis.Branch( DstarString[k], &DstarVarVec[k], "var_"+DstarString[k]+"/D"); }
#if GENvaraibles == 0	
#elif GENvaraibles == 1
	double PUWeight = 0.;
	TBranch *b_PUWeight; 
	t1->SetBranchAddress("PUWeight",&PUWeight,&b_PUWeight);
	t_analysis.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif
	//----------------------------------------
	TString htitle = "Invariant Mass of the D* ; Mass [GeV] ; Events ";
	TH1 *DsmassHisto0 = makeTH1("DsmassHisto0", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassHisto1 = makeTH1("DsmassHisto1", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassHisto2 = makeTH1("DsmassHisto2", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassHisto3 = makeTH1("DsmassHisto3", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassHisto4 = makeTH1("DsmassHisto4", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassHisto5 = makeTH1("DsmassHisto5", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassHisto6 = makeTH1("DsmassHisto6", 100,1.93,2.1, htitle, kRed);
	TH1* hDstarMass[] = { DsmassHisto0, DsmassHisto1, DsmassHisto2, DsmassHisto3, DsmassHisto4, DsmassHisto5, DsmassHisto6};
	TH1 *DsmassPtHisto1 = makeTH1("DsmassPtHisto1", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassPtHisto2 = makeTH1("DsmassPtHisto2", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassPtHisto3 = makeTH1("DsmassPtHisto3", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassPtHisto4 = makeTH1("DsmassPtHisto4", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassPtHisto5 = makeTH1("DsmassPtHisto5", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassPtHisto6 = makeTH1("DsmassPtHisto6", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassPtHisto7 = makeTH1("DsmassPtHisto7", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassPtHisto8 = makeTH1("DsmassPtHisto8", 100,1.93,2.1, htitle, kRed);
	TH1 *DsmassPtHisto9 = makeTH1("DsmassPtHisto9", 100,1.93,2.1, htitle, kRed);
	TH1* hDstarMassPt[] = { DsmassPtHisto1, DsmassPtHisto2, DsmassPtHisto3, DsmassPtHisto4, DsmassPtHisto5, DsmassPtHisto6, DsmassPtHisto7, DsmassPtHisto8, DsmassPtHisto9};

	//D* PU weight
	/*TH1 *hDsmass = makeTH1("hDsmass", 100,1.93,2.1, "Invariant Mass of the D* ; Mass [GeV/c2] ; Events ", kRed);
	TH1 *hDseta = makeTH1("hDseta", 100,-4,4, "Pseudo-rapidity distribuition of the D0(From D*) ; #eta ; Events", kRed);
	TH1 *hDsphi = makeTH1("hDsphi", 100,-4,4, "#Phi distribuition of the D* ; #Phi ; Events ", kRed);
	TH1 *hDspt = makeTH1("hDspt", 100,0,30, "pT distribuition of the D* ; pT [GeV] ; Events", kRed);
	TH1 *hD0mass = makeTH1("hD0mass", 100,1.76,1.96, "Invariant Mass of the D0(FromD*) ; Mass [GeV] ; Events ", kRed);
	TH1 *hD0eta = makeTH1("hD0eta", 100,-4,4, "Invariant Mass of the D0(FromD*) ; Mass [GeV] ; Events ", kRed);
	TH1 *hD0phi = makeTH1("hD0phi", 100,-4,4, "Invariant Mass of the D0 ; Mass [GeV] ; Events ", kRed);
	TH1 *hD0pt = makeTH1("hD0pt", 100,0,30, "Invariant Mass of the D0(FromD*) ; Mass [GeV] ; Events ", kRed);
	TH1 *hDslifetime = makeTH1("hDslifetime",100,0,2*(pow(10,-12)), "lifetime of the D0(FromD*) ; lifetime [s] ; Events ", kRed);
	TH1 *hDsMinusD0 = makeTH1("hDsMinusD0", 100,0.13,0.17, "#Delta m = m(K#pi#pi_{slow}) - m(K#pi) ; #Delta m [GeV/c2]; Events", kRed);
	//TH1* HistogramsDstar[] = {hDsmass, hDseta, hDsphi, hDspt, hD0mass, hD0eta, hD0phi, hD0pt, hDslifetime, hDsMinusD0 };*/


	//----------------------------------------
	//For D* Wrong combination
	//-----------------------------------------
	TString DstarWrongString [] = {"D0massWrong", "D0etaWrong", "D0phiWrong", "D0ptWrong", "DsmassWrong", "DsetaWrong", "DsphiWrong", "DsptWrong", };
	vector< vector<double>* > DstarVecWrong;
	vector<double>* D0massWrong =0; vector<double>* D0etaWrong =0; vector<double>* D0phiWrong =0; vector<double>* D0ptWrong =0; vector<double>* DsmassWrong =0; vector<double>* DsetaWrong =0; vector<double>* DsphiWrong =0; vector<double>* DsptWrong =0; 
	DstarVecWrong.push_back(D0massWrong); DstarVecWrong.push_back(D0etaWrong); DstarVecWrong.push_back(D0phiWrong); DstarVecWrong.push_back(D0ptWrong); DstarVecWrong.push_back(DsmassWrong); DstarVecWrong.push_back(DsetaWrong); DstarVecWrong.push_back(DsphiWrong); DstarVecWrong.push_back(DsptWrong); 
	// Input Branches
	vector< TBranch* > DstarBranchWrong; DstarBranchWrong.clear();
	TBranch *b_D0massWrong; TBranch *b_D0etaWrong; TBranch *b_D0phiWrong; TBranch *b_D0ptWrong; TBranch *b_DsmassWrong; TBranch *b_DsetaWrong; TBranch *b_DsphiWrong; TBranch *b_DsptWrong; 
	DstarBranchWrong.push_back(b_D0massWrong); DstarBranchWrong.push_back(b_D0etaWrong); DstarBranchWrong.push_back(b_D0phiWrong); DstarBranchWrong.push_back(b_D0ptWrong); DstarBranchWrong.push_back(b_DsmassWrong); DstarBranchWrong.push_back(b_DsetaWrong); DstarBranchWrong.push_back(b_DsphiWrong); DstarBranchWrong.push_back(b_DsptWrong); 
	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < DstarBranchWrong.size(); k++)
	{ t1->SetBranchAddress( DstarWrongString[k], &DstarVecWrong[k], &DstarBranchWrong[k]); }
	//-----------------------------------------
	//Variables 
	double var_D0massWrong =0.; double var_D0etaWrong =0.; double var_D0phiWrong =0.; double var_D0ptWrong =0.; double var_DsmassWrong =0.; double var_DsetaWrong =0.; double var_DsphiWrong =0.; double var_DsptWrong =0.;  
	double DstarVarWrong[] = {var_D0massWrong, var_D0etaWrong, var_D0phiWrong, var_D0ptWrong, var_DsmassWrong, var_DsetaWrong, var_DsphiWrong, var_DsptWrong, };
	sizeArray = sizeof(DstarVarWrong)/sizeof(DstarVarWrong[0]);
	vector<double> DstarVarVecWrong; DstarVarVecWrong.clear();
	DstarVarVecWrong.insert(DstarVarVecWrong.begin(), DstarVarWrong, DstarVarWrong+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < DstarVarVecWrong.size(); k++)
	{ t_DstarWrongCombination.Branch( DstarWrongString[k], &DstarVarVecWrong[k], "var_"+DstarWrongString[k]+"/D"); }


	//----------------------------------------
	//For D0
	//-----------------------------------------
	if (debug)cout << "debug D0 Variables --------------------"<< endl;
	TString D0String [] = {"D0Kpimass", "D0Kpipt", "D0Kpieta", "D0Kpiphi", "D0lifetime", "D0Kpi_VtxProb", "D0KpiDispAngle", "D0KpisXY", "D0KpidXY", "D0Kpis3D", "D0Kpid3D", "D0KpikT"};
	vector< vector<double>* > D0Vec;
	vector<double>* D0Kpimass =0; vector<double>* D0Kpipt =0; vector<double>* D0Kpieta =0; vector<double>* D0Kpiphi =0; vector<double>* D0lifetime =0; vector<double>* D0Kpi_VtxProb =0; vector<double>* D0KpiDispAngle =0; vector<double>* D0KpisXY =0; vector<double>* D0KpidXY =0; vector<double>* D0Kpis3D =0; vector<double>* D0Kpid3D =0; vector<double>* D0KpikT =0;
	D0Vec.push_back(D0Kpimass); D0Vec.push_back(D0Kpipt); D0Vec.push_back(D0Kpieta); D0Vec.push_back(D0Kpiphi); D0Vec.push_back(D0lifetime); D0Vec.push_back(D0Kpi_VtxProb); D0Vec.push_back(D0KpiDispAngle); D0Vec.push_back(D0KpisXY); D0Vec.push_back(D0KpidXY); D0Vec.push_back(D0Kpis3D); D0Vec.push_back(D0Kpid3D); D0Vec.push_back(D0KpikT);
	// Input Branches
	vector< TBranch* > D0Branch; D0Branch.clear();
	TBranch *b_D0Kpimass; TBranch *b_D0Kpipt; TBranch *b_D0Kpieta; TBranch *b_D0Kpiphi; TBranch *b_D0lifetime; TBranch *b_D0Kpi_VtxProb; TBranch *b_D0KpiDispAngle; TBranch *b_D0KpisXY; TBranch *b_D0KpidXY; TBranch *b_D0Kpis3D; TBranch *b_D0Kpid3D; TBranch *b_D0KpikT; 
	D0Branch.push_back(b_D0Kpimass); D0Branch.push_back(b_D0Kpipt); D0Branch.push_back(b_D0Kpieta); D0Branch.push_back(b_D0Kpiphi); D0Branch.push_back(b_D0lifetime); D0Branch.push_back(b_D0Kpi_VtxProb); D0Branch.push_back(b_D0KpiDispAngle); D0Branch.push_back(b_D0KpisXY); D0Branch.push_back(b_D0KpidXY); D0Branch.push_back(b_D0Kpis3D); D0Branch.push_back(b_D0Kpid3D); D0Branch.push_back(b_D0KpikT);
	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < D0Branch.size(); k++)
	{ t1->SetBranchAddress( D0String[k], &D0Vec[k], &D0Branch[k]); }
	//-----------------------------------------
	//Variables
	double var_D0Kpimass =0.; double var_D0Kpipt =0.; double var_D0Kpieta =0.; double var_D0Kpiphi =0.; double var_D0lifetime =0.; double var_D0Kpi_VtxProb =0.; double var_D0KpiDispAngle =0.; double var_D0KpisXY =0.; double var_D0KpidXY =0.; double var_D0Kpis3D =0.; double var_D0Kpid3D =0.; double var_D0KpikT =0.;  
	double D0Var[] = {var_D0Kpimass, var_D0Kpipt, var_D0Kpieta, var_D0Kpiphi, var_D0lifetime, var_D0Kpi_VtxProb, var_D0KpiDispAngle, var_D0KpisXY, var_D0KpidXY, var_D0Kpis3D, var_D0Kpid3D, var_D0KpikT };
	sizeArray = sizeof(D0Var)/sizeof(D0Var[0]); // Number of elements
	vector<double> D0VarVec; D0VarVec.clear();
	D0VarVec.insert(D0VarVec.begin(), D0Var, D0Var+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < D0VarVec.size(); k++)
	{ t_D0analysis.Branch( D0String[k], &D0VarVec[k], "var_"+D0String[k]+"/D"); }
#if GENvaraibles == 0	
#elif GENvaraibles == 1
	t_D0analysis.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif

	//Study of Mass with different D0 significance cuts
	htitle = "Invariant Mass of the D0 ; Mass [GeV] ; Events ";
	TH1 *D0KpimassHisto0 = makeTH1("D0KpimassHisto0", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto1 = makeTH1("D0KpimassHisto1", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto2 = makeTH1("D0KpimassHisto2", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto3 = makeTH1("D0KpimassHisto3", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto4 = makeTH1("D0KpimassHisto4", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto5 = makeTH1("D0KpimassHisto5", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto6 = makeTH1("D0KpimassHisto6", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto7 = makeTH1("D0KpimassHisto7", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto8 = makeTH1("D0KpimassHisto8", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto9 = makeTH1("D0KpimassHisto9", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassHisto10 = makeTH1("D0KpimassHisto10", 100,1.76,1.96, htitle, kRed);
	TH1* hD0Mass[] = { D0KpimassHisto0, D0KpimassHisto1, D0KpimassHisto2, D0KpimassHisto3, D0KpimassHisto4, D0KpimassHisto5, D0KpimassHisto6, D0KpimassHisto7, D0KpimassHisto8, D0KpimassHisto9, D0KpimassHisto10};

	TH1 *D0KpimassPtHisto1 = makeTH1("D0KpimassPtHisto1", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassPtHisto2 = makeTH1("D0KpimassPtHisto2", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassPtHisto3 = makeTH1("D0KpimassPtHisto3", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassPtHisto4 = makeTH1("D0KpimassPtHisto4", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassPtHisto5 = makeTH1("D0KpimassPtHisto5", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassPtHisto6 = makeTH1("D0KpimassPtHisto6", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassPtHisto7 = makeTH1("D0KpimassPtHisto7", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassPtHisto8 = makeTH1("D0KpimassPtHisto8", 100,1.76,1.96, htitle, kRed);
	TH1 *D0KpimassPtHisto9 = makeTH1("D0KpimassPtHisto9", 100,1.76,1.96, htitle, kRed);
	TH1* hD0MassPt[] = { D0KpimassPtHisto1, D0KpimassPtHisto2, D0KpimassPtHisto3, D0KpimassPtHisto4, D0KpimassPtHisto5, D0KpimassPtHisto6, D0KpimassPtHisto7, D0KpimassPtHisto8, D0KpimassPtHisto9 };

#if GENvaraibles == 0	
#elif GENvaraibles == 1
	//----------------------------------------
	//For D* MC
	//-----------------------------------------
	if (debug)cout << "debug D* MC Variables --------------------"<< endl;
	TString MCDstarString [] = {"MCDseta", "MCDsphi", "MCDspt", "MCDsenergy", "MCDsp", "MCDset", "MCDsmass", "MCD0eta", "MCD0phi", "MCD0pt", "MCD0energy", "MCD0p", "MCD0et", "MCD0rapidity", "MCD0mass", "MCD0lifetime" };
	vector< vector<double>* > MCDstarVec;
	vector<double>* MCDseta =0; vector<double>* MCDsphi =0; vector<double>* MCDspt =0; vector<double>* MCDsenergy =0; vector<double>* MCDsp =0; vector<double>* MCDset =0; vector<double>* MCDsmass =0; vector<double>* MCD0eta =0; vector<double>* MCD0phi =0; vector<double>* MCD0pt =0; vector<double>* MCD0energy =0; vector<double>* MCD0p =0; vector<double>* MCD0et =0; vector<double>* MCD0rapidity =0; vector<double>* MCD0mass =0; vector<double>* MCD0lifetime =0;
	MCDstarVec.push_back(MCDseta); MCDstarVec.push_back(MCDsphi); MCDstarVec.push_back(MCDspt); MCDstarVec.push_back(MCDsenergy); MCDstarVec.push_back(MCDsp); MCDstarVec.push_back(MCDset); MCDstarVec.push_back(MCDsmass); MCDstarVec.push_back(MCD0eta); MCDstarVec.push_back(MCD0phi); MCDstarVec.push_back(MCD0pt); MCDstarVec.push_back(MCD0energy); MCDstarVec.push_back(MCD0p); MCDstarVec.push_back(MCD0et); MCDstarVec.push_back(MCD0rapidity); MCDstarVec.push_back(MCD0mass); MCDstarVec.push_back(MCD0lifetime);
	// Input Branches
	vector< TBranch* > MCDstarBranch; MCDstarBranch.clear();
	TBranch *b_MCDseta; TBranch *b_MCDsphi; TBranch *b_MCDspt; TBranch *b_MCDsenergy; TBranch *b_MCDsp; TBranch *b_MCDset; TBranch *b_MCDsmass; TBranch *b_MCD0eta; TBranch *b_MCD0phi; TBranch *b_MCD0pt; TBranch *b_MCD0energy; TBranch *b_MCD0p; TBranch *b_MCD0et; TBranch *b_MCD0rapidity; TBranch *b_MCD0mass; TBranch *b_MCD0lifetime; 
	MCDstarBranch.push_back(b_MCDseta); MCDstarBranch.push_back(b_MCDsphi); MCDstarBranch.push_back(b_MCDspt); MCDstarBranch.push_back(b_MCDsenergy); MCDstarBranch.push_back(b_MCDsp); MCDstarBranch.push_back(b_MCDset); MCDstarBranch.push_back(b_MCDsmass); MCDstarBranch.push_back(b_MCD0eta); MCDstarBranch.push_back(b_MCD0phi); MCDstarBranch.push_back(b_MCD0pt); MCDstarBranch.push_back(b_MCD0energy); MCDstarBranch.push_back(b_MCD0p); MCDstarBranch.push_back(b_MCD0et); MCDstarBranch.push_back(b_MCD0rapidity); MCDstarBranch.push_back(b_MCD0mass); MCDstarBranch.push_back(b_MCD0lifetime);
	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < MCDstarBranch.size(); k++)
	{ t1->SetBranchAddress( MCDstarString[k], &MCDstarVec[k], &MCDstarBranch[k]); }
	//-----------------------------------------
	//Variables 
	double var_MCDseta =0.; double var_MCDsphi =0.; double var_MCDspt =0.; double var_MCDsenergy =0.; double var_MCDsp =0.; double var_MCDset =0.; double var_MCDsmass =0.; double var_MCD0eta =0.; double var_MCD0phi =0.; double var_MCD0pt =0.; double var_MCD0energy =0.; double var_MCD0p =0.; double var_MCD0et =0.; double var_MCD0rapidity =0.; double var_MCD0mass =0.; double var_MCD0lifetime =0.;  
	double MCDstarVar[] = {var_MCDseta, var_MCDsphi, var_MCDspt, var_MCDsenergy, var_MCDsp, var_MCDset, var_MCDsmass, var_MCD0eta, var_MCD0phi, var_MCD0pt, var_MCD0energy, var_MCD0p, var_MCD0et, var_MCD0rapidity, var_MCD0mass, var_MCD0lifetime };
	sizeArray = sizeof(MCDstarVar)/sizeof(MCDstarVar[0]); // Number of elements
	vector<double> MCDstarVarVec; MCDstarVarVec.clear();
	MCDstarVarVec.insert(MCDstarVarVec.begin(), MCDstarVar, MCDstarVar+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < MCDstarVarVec.size(); k++)
	{ t_DsMC.Branch( MCDstarString[k], &MCDstarVarVec[k], "var_"+MCDstarString[k]+"/D"); }
	t_DsMC.Branch("PUWeight",&PUWeight, "PUWeight/D");
	//----------------------------------------
	//For D* Contamination
	vector<int>* FlagDstarfromB = 0;
	TBranch *b_FlagDstarfromB; t1->SetBranchAddress("FlagDstarfromB",&FlagDstarfromB,&b_FlagDstarfromB);
	double var_DstarFromBmass = 0.; double var_DstarFromBpt = 0.; double var_DstarFromBeta = 0.;
	double var_DstarFromPPmass = 0.;	double var_DstarFromPPpt = 0.; double var_DstarFromPPeta = 0.;
	int DstarSec =0; int DstarPrompt =0;
	t_DstarContamination.Branch("DstarFromBmass",&var_DstarFromBmass, "var_DstarFromBmass/D");
	t_DstarContamination.Branch("DstarFromBpt",&var_DstarFromBpt, "var_DstarFromBpt/D");
	t_DstarContamination.Branch("DstarFromBeta",&var_DstarFromBeta, "var_DstarFromBeta/D");
	t_DstarContamination.Branch("DstarFromPPmass",&var_DstarFromPPmass, "var_DstarFromPPmass/D");
	t_DstarContamination.Branch("DstarFromPPpt",&var_DstarFromPPpt, "var_DstarFromPPpt/D");
	t_DstarContamination.Branch("DstarFromPPeta",&var_DstarFromPPeta, "var_DstarFromPPeta/D");
	
	//----------------------------------------
	// D* reconstruction Efficiency 
	unsigned long GenDspT1, GenDspT2, GenDspT3, GenDspT4, GenDspT5, GenDspT6, GenDspT7, GenDspT8, GenDspT9;
	GenDspT1 = 0; GenDspT2 = 0; GenDspT3 = 0; GenDspT4 = 0; GenDspT5 = 0; GenDspT6 = 0; GenDspT7 = 0; GenDspT8 = 0; GenDspT9 = 0;
	unsigned long GenDspT[] = { GenDspT1, GenDspT2, GenDspT3, GenDspT4, GenDspT5, GenDspT6, GenDspT7, GenDspT8, GenDspT9 };

	unsigned long RecDspT1, RecDspT2, RecDspT3, RecDspT4, RecDspT5, RecDspT6, RecDspT7, RecDspT8, RecDspT9;
	RecDspT1 = 0; RecDspT2 = 0; RecDspT3 = 0; RecDspT4 = 0; RecDspT5 = 0; RecDspT6 = 0; RecDspT7 = 0; RecDspT8 = 0; RecDspT9 = 0;
	unsigned long RecDspT[] = { RecDspT1, RecDspT2, RecDspT3, RecDspT4, RecDspT5, RecDspT6, RecDspT7, RecDspT8, RecDspT9};
	
	unsigned long GenDsEta1, GenDsEta2, GenDsEta3, GenDsEta4, GenDsEta5, GenDsEta6, GenDsEta7, GenDsEta8, GenDsEta9, GenDsEta10;
	GenDsEta1 = 0; GenDsEta2 = 0; GenDsEta3 = 0; GenDsEta4 = 0; GenDsEta5 = 0; GenDsEta6 = 0; GenDsEta7 = 0; GenDsEta8 = 0; GenDsEta9 = 0; GenDsEta10 = 0;
	unsigned long GenDsEta[] = { GenDsEta1, GenDsEta2, GenDsEta3, GenDsEta4, GenDsEta5, GenDsEta6, GenDsEta7, GenDsEta8, GenDsEta9, GenDsEta10 };

	unsigned long RecDsEta1, RecDsEta2, RecDsEta3, RecDsEta4, RecDsEta5, RecDsEta6, RecDsEta7, RecDsEta8, RecDsEta9, RecDsEta10;
	RecDsEta1 = 0; RecDsEta2 = 0; RecDsEta3 = 0; RecDsEta4 = 0; RecDsEta5 = 0; RecDsEta6 = 0; RecDsEta7 = 0; RecDsEta8 = 0; RecDsEta9 = 0; RecDsEta10 = 0;
	unsigned long RecDsEta[] = { RecDsEta1, RecDsEta2, RecDsEta3, RecDsEta4, RecDsEta5, RecDsEta6, RecDsEta7, RecDsEta8, RecDsEta9, RecDsEta10 };

	TH2* hdeltaRDstarMass = makeTH2GenRec( "hdeltaRDstarMass", "Mass[GeV]", 100, 1.94 , 2.1 );
	TH2* hdeltaRDstarPt = makeTH2GenRec( "hdeltaRDstarPt", "pT[GeV]", 100, -1, 40. );
	TH2* hdeltaRDstarEta = makeTH2GenRec( "hdeltaRDstarEta", "#eta", 100, -1, 3. );
	TH2* hdeltaRDstarTime = makeTH2GenRec( "hdeltaRDstarTime", "t[s]", 100, 0, 1*pow(10,-12) );

	TH2* hdeltaRD0Mass = makeTH2GenRec( "hdeltaRD0Mass", "Mass[GeV]", 100, 1.76, 1.96 );
	TH2* hdeltaRD0Pt = makeTH2GenRec( "hdeltaRD0Pt", "pT[GeV]", 100, -1, 30. );
	TH2* hdeltaRD0Eta = makeTH2GenRec( "hdeltaRD0Eta", "#eta", 100, -1, 3. );
	TH2* hdeltaRD0Time = makeTH2GenRec( "hdeltaRD0Time", "t[s]", 100, 0, 3*pow(10,-12) );

	
	const int nbinsPt_left = 9;
	double binsPt_left[nbinsPt_left+1] = {4., 5., 6., 7., 8., 12., 16., 24., 40., 100.};
	const int nbinsEta_left = 10;
	double binsEta_left[nbinsEta_left+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1};
	TH1 *hRecEta_MC = makeTH1Rec("hRecEta_MC", nbinsEta_left, binsEta_left );
	TH1 *hRecPt_MC = makeTH1Rec("hRecPt_MC", nbinsPt_left, binsPt_left );
	TH1 *hRecEta_Rec1 = makeTH1Rec("hRecEta_Rec1", nbinsEta_left, binsEta_left );
	TH1 *hRecPt_Rec1 = makeTH1Rec("hRecPt_Rec1", nbinsPt_left, binsPt_left );


	
	//----------------------------------------
	//For D0 MC
	//-----------------------------------------
	TString MCD0String [] = {"MCpromptD0eta", "MCpromptD0phi", "MCpromptD0pt", "MCpromptD0energy", "MCpromptD0p", "MCpromptD0et", "MCpromptD0rapidity", "MCpromptD0mass", "MCpromptD0_DispAngle", "MCpromptD0lifetime" };
	vector< vector<double>* > MCD0Vec;
	vector<double>* MCpromptD0eta =0; vector<double>* MCpromptD0phi =0; vector<double>* MCpromptD0pt =0; vector<double>* MCpromptD0energy =0; vector<double>* MCpromptD0p =0; vector<double>* MCpromptD0et =0; vector<double>* MCpromptD0rapidity =0; vector<double>* MCpromptD0mass =0; vector<double>* MCpromptD0_DispAngle =0; vector<double>* MCpromptD0lifetime =0; 
	MCD0Vec.push_back(MCpromptD0eta); MCD0Vec.push_back(MCpromptD0phi); MCD0Vec.push_back(MCpromptD0pt); MCD0Vec.push_back(MCpromptD0energy); MCD0Vec.push_back(MCpromptD0p); MCD0Vec.push_back(MCpromptD0et); MCD0Vec.push_back(MCpromptD0rapidity); MCD0Vec.push_back(MCpromptD0mass); MCD0Vec.push_back(MCpromptD0_DispAngle); MCD0Vec.push_back(MCpromptD0lifetime); 
	// Input Branches
	vector< TBranch* > MCD0Branch; MCD0Branch.clear();
	TBranch *b_MCpromptD0eta; TBranch *b_MCpromptD0phi; TBranch *b_MCpromptD0pt; TBranch *b_MCpromptD0energy; TBranch *b_MCpromptD0p; TBranch *b_MCpromptD0et; TBranch *b_MCpromptD0rapidity; TBranch *b_MCpromptD0mass; TBranch *b_MCpromptD0_DispAngle; TBranch *b_MCpromptD0lifetime; 
	MCD0Branch.push_back(b_MCpromptD0eta); MCD0Branch.push_back(b_MCpromptD0phi); MCD0Branch.push_back(b_MCpromptD0pt); MCD0Branch.push_back(b_MCpromptD0energy); MCD0Branch.push_back(b_MCpromptD0p); MCD0Branch.push_back(b_MCpromptD0et); MCD0Branch.push_back(b_MCpromptD0rapidity); MCD0Branch.push_back(b_MCpromptD0mass); MCD0Branch.push_back(b_MCpromptD0_DispAngle); MCD0Branch.push_back(b_MCpromptD0lifetime); 
	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < MCD0Branch.size(); k++)
	{ t1->SetBranchAddress( MCD0String[k], &MCD0Vec[k], &MCD0Branch[k]); }
	//-----------------------------------------
	//Variables 
	double var_MCpromptD0eta =0.; double var_MCpromptD0phi =0.; double var_MCpromptD0pt =0.; double var_MCpromptD0energy =0.; double var_MCpromptD0p =0.; double var_MCpromptD0et =0.; double var_MCpromptD0rapidity =0.; double var_MCpromptD0mass =0.; double var_MCpromptD0_DispAngle =0.; double var_MCpromptD0lifetime =0.;  
	double MCD0Var[] = {var_MCpromptD0eta, var_MCpromptD0phi, var_MCpromptD0pt, var_MCpromptD0energy, var_MCpromptD0p, var_MCpromptD0et, var_MCpromptD0rapidity, var_MCpromptD0mass, var_MCpromptD0_DispAngle, var_MCpromptD0lifetime };
	sizeArray = sizeof(MCD0Var)/sizeof(MCD0Var[0]); // Number of elements
	vector<double> MCD0VarVec; MCD0VarVec.clear();
	MCD0VarVec.insert(MCD0VarVec.begin(), MCD0Var, MCD0Var+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < MCD0VarVec.size(); k++)
	{ t_D0MC.Branch( MCD0String[k], &MCD0VarVec[k], "var_"+MCD0String[k]+"/D"); }
	t_D0MC.Branch("PUWeight",&PUWeight, "PUWeight/D");
	//----------------------------------------
	//For D0 Contamination
	vector<int>* FlagD0fromB = 0;
	TBranch *b_FlagD0fromB; t1->SetBranchAddress("FlagD0fromB",&FlagD0fromB,&b_FlagD0fromB);
	double var_D0FromBmass = 0.; double var_D0FromBpt = 0.; double var_D0FromBeta = 0.;
	double var_D0FromPPmass = 0.;	double var_D0FromPPpt = 0.; double var_D0FromPPeta = 0.;
	int D0Sec =0; int D0Prompt =0;
	t_D0Contamination.Branch("D0FromBmass",&var_D0FromBmass, "var_D0FromBmass/D");
	t_D0Contamination.Branch("D0FromBpt",&var_D0FromBpt, "var_D0FromBpt/D");
	t_D0Contamination.Branch("D0FromBeta",&var_D0FromBeta, "var_D0FromBeta/D");
	t_D0Contamination.Branch("D0FromPPmass",&var_D0FromPPmass, "var_D0FromPPmass/D");
	t_D0Contamination.Branch("D0FromPPpt",&var_D0FromPPpt, "var_D0FromPPpt/D");
	t_D0Contamination.Branch("D0FromPPeta",&var_D0FromPPeta, "var_D0FromPPeta/D");
	//----------------------------------------
	// D0 recosntruction Efficiency 
	unsigned long GenD0pT1, GenD0pT2, GenD0pT3, GenD0pT4, GenD0pT5, GenD0pT6, GenD0pT7, GenD0pT8, GenD0pT9;
	GenD0pT1 = 0; GenD0pT2 = 0; GenD0pT3 = 0; GenD0pT4 = 0; GenD0pT5 = 0; GenD0pT6 = 0; GenD0pT7 = 0; GenD0pT8 = 0; GenD0pT9 = 0;
	unsigned long GenD0pT[] = { GenD0pT1, GenD0pT2, GenD0pT3, GenD0pT4, GenD0pT5, GenD0pT6, GenD0pT7, GenD0pT8, GenD0pT9 };

	unsigned long RecD0pT1, RecD0pT2, RecD0pT3, RecD0pT4, RecD0pT5, RecD0pT6, RecD0pT7, RecD0pT8, RecD0pT9;
	RecD0pT1 = 0; RecD0pT2 = 0; RecD0pT3 = 0; RecD0pT4 = 0; RecD0pT5 = 0; RecD0pT6 = 0; RecD0pT7 = 0; RecD0pT8 = 0; RecD0pT9 = 0;
	unsigned long RecD0pT[] = { RecD0pT1, RecD0pT2, RecD0pT3, RecD0pT4, RecD0pT5, RecD0pT6, RecD0pT7, RecD0pT8, RecD0pT9};
	
	unsigned long GenD0Eta1, GenD0Eta2, GenD0Eta3, GenD0Eta4, GenD0Eta5, GenD0Eta6, GenD0Eta7, GenD0Eta8, GenD0Eta9, GenD0Eta10;
	GenD0Eta1 = 0; GenD0Eta2 = 0; GenD0Eta3 = 0; GenD0Eta4 = 0; GenD0Eta5 = 0; GenD0Eta6 = 0; GenD0Eta7 = 0; GenD0Eta8 = 0; GenD0Eta9 = 0; GenD0Eta10 = 0;
	unsigned long GenD0Eta[] = { GenD0Eta1, GenD0Eta2, GenD0Eta3, GenD0Eta4, GenD0Eta5, GenD0Eta6, GenD0Eta7, GenD0Eta8, GenD0Eta9, GenD0Eta10 };

	unsigned long RecD0Eta1, RecD0Eta2, RecD0Eta3, RecD0Eta4, RecD0Eta5, RecD0Eta6, RecD0Eta7, RecD0Eta8, RecD0Eta9, RecD0Eta10;
	RecD0Eta1 = 0; RecD0Eta2 = 0; RecD0Eta3 = 0; RecD0Eta4 = 0; RecD0Eta5 = 0; RecD0Eta6 = 0; RecD0Eta7 = 0; RecD0Eta8 = 0; RecD0Eta9 = 0; RecD0Eta10 = 0;
	unsigned long RecD0Eta[] = { RecD0Eta1, RecD0Eta2, RecD0Eta3, RecD0Eta4, RecD0Eta5, RecD0Eta6, RecD0Eta7, RecD0Eta8, RecD0Eta9, RecD0Eta10 };
#endif
	//----------------------------------------
	// LOOP TREE ENTRIES FOR FILE f1
	for (Long64_t jentry=0; jentry < nentries; jentry++) 
	{
		if (debug)cout << "debug File loop --------------------"<< endl;
		Long64_t ientry = t1->LoadTree(jentry);
      if (ientry < 0) break;
		//-----------------------------------------
		//Output about percent program executed
		double percent = (jentry*100)/nentries;		
		if ( jentry % partpercent == 0) cout <<"*****"<< percent << "per cent done***" << endl;
		//-----------------------------------------
		// Get the i entry
		//----------------------------------------
		//For D*
		for( unsigned int k = 0; k < DstarVarVec.size(); k++)
		{ DstarBranch[k]->GetEntry(ientry); }
		//----------------------------------------
		//For D* wrong combination
		for( unsigned int k = 0; k < DstarVarVecWrong.size(); k++)
		{ DstarBranchWrong[k]->GetEntry(ientry); }
		//----------------------------------------
		//For D0
		for( unsigned int k = 0; k < D0VarVec.size(); k++)
		{ D0Branch[k]->GetEntry(ientry); }
#if GENvaraibles == 0	
#elif GENvaraibles == 1
		b_PUWeight->GetEntry(ientry);
		//cout << "PUWeight: " << PUWeight << endl;
		//----------------------------------------
		//For D* MC
		for( unsigned int k = 0; k < MCDstarVarVec.size(); k++)
		{ MCDstarBranch[k]->GetEntry(ientry); }
		b_FlagDstarfromB->GetEntry(ientry);
		//----------------------------------------
		//For D* MC
		for( unsigned int k = 0; k < MCD0VarVec.size(); k++)
		{ MCD0Branch[k]->GetEntry(ientry); }
		b_FlagD0fromB->GetEntry(ientry);
#endif
		
		//----------------------------------------
		// SELECTION
		//----------------------------------------
		//For D* and D0(from D*)
		if ( (DstarVec[0]->size() > 0) && FlagDs){
		for(unsigned int i=0; i < DstarVec[0]->size(); i++)
		{  
			//{"Dsmass", "Dseta", "Dsphi", "Dspt", "DSDeltaR", "D0mass", "D0eta", "D0phi", "D0pt", "D0fromDSdXY", "D0fromDSsXY", "D0fromDSs3D", "D0fromDSd3D", "Anglephi", "D0_VtxProb"};
			//PossibleDs++;					
			if( DstarVec[14]->at(i) < 0.01) continue; //DsD0ProbVtx++;		
			if( DstarVec[13]->at(i) < 0.99 ) continue; //DsD0Angle++;	
			if( fabs(DstarVec[5]->at(i) - 1.86484) > 0.1 ) continue; //D0fromDsMinusPDG++;
			if( DstarVec[8]->at(i) < 3. ) continue; //DsD0pt++;
			if( (DstarVec[0]->at(i) - DstarVec[5]->at(i)) > 0.16) continue; //DsD0deltaM++;*/
		
			//Evolution Significance Cut
			double binDstarMass[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
			sizeArray = sizeof(hDstarMass)/sizeof(hDstarMass[0]); // Number of elements
			for( int k = 0; k < sizeArray; k++)
			{ if( DstarVec[11]->at(i) > binDstarMass[k]) { hDstarMass[k]->Fill(DstarVec[0]->at(i));}	}

			if( DstarVec[11]->at(i) < SigCutD0fromDs) continue; //DsD0Sig++;

			//Evolution pT Cut
			double binDstarMassPt[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
			sizeArray = sizeof(hDstarMassPt)/sizeof(hDstarMassPt[0]); // Number of elements
			for( int k = 0; k < sizeArray; k++)
			{ if( DstarVec[3]->at(i) > binDstarMassPt[k]) { hDstarMassPt[k]->Fill(DstarVec[0]->at(i));}	}

			//----------------------------------------------------
			// Fill Branches
			for( unsigned int k = 0; k < DstarVarVec.size(); k++)
			{	DstarVarVec[k] = DstarVec[k]->at(i); }

#if GENvaraibles == 0
			// Fill histograms
			//for( unsigned int k = 0; k < DstarVarVec.size(); k++)
			//{ HistogramsDstar[k]->Fill(DstarVec[0]->at(i));	}
#elif GENvaraibles == 1
			// Fill histograms
			//for( unsigned int k = 0; k < DstarVarVec.size(); k++)
			//{ HistogramsDstar[k]->Fill( DstarVec[0]->at(i), PUWeight);	}
#endif

			//TH1* HistogramsDstar[] = {hDsmass, hDseta, hDsphi, hDspt, hD0mass, hD0eta, hD0phi, hD0pt, hDslifetime, hDsMinusD0 };
			//hDslifetime->Fill(( D0mass->at(i)*D0fromDSdXY->at(i)*pow(10,-2))/(c*D0pt->at(i)));
			//hDsMinusD0->Fill(( Dsmass->at(i) - D0mass->at(i)));



			//cout << "Dsmass: "<<  DstarVec[0]->at(i) << " # D0mass: "<<  DstarVec[5]->at(i) << " # Prob: "<<  DstarVec[14]->at(i) << " # AnglePhi: "<<  DstarVec[13]->at(i) << endl;	
			t_analysis.Fill();

			if (debug)cout << "debug D* 12 --------------------"<< endl;
	  	}}
		//----------------------------------------
		//For D0 
		//-----------------------------------------
		if ( (D0Vec[0]->size() > 0) && FlagD0){
		for(unsigned int i=0; i < D0Vec[0]->size(); i++)
		{	

			//{"D0Kpimass", "D0Kpipt", "D0Kpieta", "D0Kpiphi", "D0lifetime", "D0Kpi_VtxProb", "D0KpiDispAngle", "D0KpisXY", "D0KpidXY", "D0Kpis3D", "D0Kpid3D", "D0KpikT"}			
			//PossibleD0++;
			if( D0Vec[5]->at(i) < 0.01) continue; //D0ProbVtx++;	
			if( D0Vec[6]->at(i) < 0.99 ) continue; //D0Angle++;
			if( abs(D0Vec[0]->at(i)-1.86484) > 0.1) continue;	//CountD0minusPDG++;

			//Evolution pT Cut
			if( D0Vec[9]->at(i) < SigCutD0)
			{	double binD0Pt[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
				sizeArray = sizeof(hDstarMassPt)/sizeof(hDstarMassPt[0]); // Number of elements
				for( int k = 0; k < sizeArray; k++)
				{ if( D0Vec[1]->at(i)> binD0Pt[k]) { hD0MassPt[k]->Fill(D0Vec[0]->at(i));}	}	
			}

			if( D0Vec[1]->at(i) < 2. ) continue; //CountD0pt++;

			//Evolution Significance Cut
			double binD0Sig[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
			sizeArray = sizeof(hDstarMassPt)/sizeof(hDstarMassPt[0]); // Number of elements
			for( int k = 0; k < sizeArray; k++)
			{ if( D0Vec[9]->at(i)>binD0Sig[k]) { hD0Mass[k]->Fill(D0Vec[0]->at(i));}	}	
			
			if( D0Vec[9]->at(i) < SigCutD0) continue; //D0Sig++;

			//----------------------------------------------------
			// Fill Branches
			for( unsigned int k = 0; k < D0VarVec.size(); k++)
			{ D0VarVec[k] = D0Vec[k]->at(i); }
		
			t_D0analysis.Fill();
		}}

#if GENvaraibles == 0	
#elif GENvaraibles == 1
		//----------------------------------------
		//For D* MC	
		if ( (MCDstarVec[0]->size() > 0) && FlagDsMC){ //D* MC protection
		for(unsigned int i=0; i < MCDstarVec[0]->size(); i++) // loop D* MC
		{
			// Fill Branches
			for( unsigned int k = 0; k < MCDstarVarVec.size(); k++)
			{ MCDstarVarVec[k] = MCDstarVec[k]->at(i); }

			hRecEta_MC->Fill(abs(MCDseta->at(i)));
			hRecPt_MC->Fill(MCDspt->at(i));

			// D* recosntruction Efficiency pT
			double DstarRangePT[] = {4., 5., 6., 7., 8., 12., 16., 24., 40., 100.};
			sizeArray = sizeof(DstarRangePT)/sizeof(DstarRangePT[0]); // Number of elements
			for( int k = 0; k < (sizeArray-1) ; k++)
			{	if( (MCDstarVec[2]->at(i) > DstarRangePT[k]) and (MCDstarVec[2]->at(i) < DstarRangePT[k+1]) )
				{	GenDspT[k]++;} 
			}

			// D* recosntruction Efficiency Eta
			double DstarRangeEta[] = {0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.1};
			sizeArray = sizeof(DstarRangeEta)/sizeof(DstarRangeEta[0]); // Number of elements
			for( int k = 0; k < (sizeArray-1) ; k++)
			{	if( (MCDstarVec[0]->at(i) > DstarRangeEta[k]) and (MCDstarVec[0]->at(i) < DstarRangeEta[k+1]) )
				{	GenDsEta[k]++;} 
			}

			if (debug)cout << "debug D* 18 --------------------"<< endl;
			
			t_DsMC.Fill();			
			//----------------------------------------
			//For D* Matching
			//----------------------------------------
			if ( DstarVec[0]->size() > 0 ){

			double minDeltaR = 1000.;
			int id = -1; // Protection
			double deltaR = 0.;
			double deltaEta= 0.;
			double deltaPhi= 0.;
			for(unsigned int j=0; j < DstarVec[0]->size(); j++)
			{	
			//{"Dsmass", "Dseta", "Dsphi", "Dspt", "DSDeltaR", "D0mass", "D0eta", "D0phi", "D0pt", "D0fromDSdXY", "D0fromDSsXY", "D0fromDSs3D", "D0fromDSd3D", "Anglephi", "D0_VtxProb"};
				if( DstarVec[14]->at(j) < 0.01) continue;	
				if( DstarVec[13]->at(j) < 0.99 ) continue; 
				if( fabs(DstarVec[5]->at(j) - 1.86484) > 0.1 ) continue;
				if( DstarVec[8]->at(j) < 3. ) continue;
				if( (DstarVec[0]->at(j) - DstarVec[5]->at(j)) > 0.16) continue;
				if( DstarVec[11]->at(j) < SigCutD0fromDs) continue; //DsD0Sig++;
						
				//----------------------------------------
				//For D* Contamination			
				if (debug)cout << "For D* Contamination"<< endl;

				if( FlagDstarfromB->size() > 0)
				{	if( FlagDstarfromB->at(i) == 1)
					{	DstarSec++;
						var_DstarFromBmass = DstarVec[0]->at(j);
						var_DstarFromBpt = DstarVec[3]->at(j);
						var_DstarFromBeta = DstarVec[1]->at(j);
					}
					else
					{	DstarPrompt++;
						var_DstarFromPPmass = DstarVec[0]->at(j);
						var_DstarFromPPpt = DstarVec[3]->at(j);
						var_DstarFromPPeta = DstarVec[1]->at(j);
					}
				}
				t_DstarContamination.Fill();

				//-----------------------------------
				//DeltaR Function
				deltaEta = MCDstarVec[0]->at(i) - DstarVec[1]->at(j);
				deltaPhi = MCDstarVec[1]->at(i) - DstarVec[2]->at(j);
				deltaR = sqrt( pow( deltaEta, 2) + pow( deltaPhi ,2));
				if( deltaR < minDeltaR) { id = j; minDeltaR = deltaR;} //Get the lower delta
				if (debug)cout << "debug D* 20 --------------------"<< endl;
			}

			if (id != -1){
				double recMass = -1; double GenMass = -1; // Protection
				double recPt = -1; double GenPt = -1;
				double recEta = -1; double GenEta = -1;
				double recTime = -1; double GenTime = -1;
				//Histograms for efficiency
				if( deltaR < 0.03){

					hRecEta_Rec1->Fill(abs(Dseta->at(id)));
					hRecPt_Rec1->Fill(Dspt->at(id));
			
					recMass = DstarVec[0]->at(id); GenMass = MCDstarVec[6]->at(i);
					recPt = DstarVec[3]->at(id); GenPt = MCDstarVec[2]->at(i);
					recEta = DstarVec[1]->at(id); GenEta = MCDstarVec[0]->at(i);
					recTime = ( (DstarVec[0]->at(id))*(DstarVec[9]->at(id))*pow(10,-2) ) / (c*(DstarVec[3]->at(id)) ); 
					GenTime = MCDstarVec[15]->at(i);

					hdeltaRDstarMass->Fill( recMass, GenMass );
					hdeltaRDstarPt->Fill( recPt, GenPt );
					hdeltaRDstarEta->Fill( abs(recEta), abs(GenEta) );
					hdeltaRDstarTime->Fill( recTime, GenTime );

					// D* recosntruction Efficiency pT
					sizeArray = sizeof(DstarRangePT)/sizeof(DstarRangePT[0]); // Number of elements
					for( int k = 0; k < (sizeArray-1) ; k++)
					{	if( (DstarVec[3]->at(id) > DstarRangePT[k] ) and ( DstarVec[3]->at(id) < DstarRangePT[k+1]) )
						{	RecDspT[k]++;} }
					// D* recosntruction Efficiency Eta
					sizeArray = sizeof(DstarRangeEta)/sizeof(DstarRangeEta[0]); // Number of elements
					for( int k = 0; k < (sizeArray-1) ; k++)
					{	if( (DstarVec[1]->at(id) > DstarRangeEta[k] ) and (DstarVec[1]->at(id) < DstarRangeEta[k+1]) )
						{	RecDsEta[k]++;} }		
				}
		
			}
			} // End D* Matching

		} // End loop D* MC
		} // End D* MC protection
		//----------------------------------------
		//For D0 MC
		if ( (MCD0Vec[0]->size() > 0) && FlagD0MC){ //D0 MC protection
		for(unsigned int i=0; i < MCD0Vec[0]->size(); i++) // loop D0 MC
		{
			// Fill Branches
			for( unsigned int k = 0; k < MCD0VarVec.size(); k++)
			{ MCD0VarVec[k] = MCD0Vec[k]->at(i); }

			// D0 recosntruction Efficiency pT
			double D0RangePT[] = {4., 5., 6., 7., 8., 12., 16., 24., 40., 100.};
			sizeArray = sizeof(D0RangePT)/sizeof(D0RangePT[0]); // Number of elements
			for( int k = 0; k < (sizeArray-1) ; k++)
			{	if( (MCD0Vec[2]->at(i) > D0RangePT[k]) and (MCD0Vec[2]->at(i) < D0RangePT[k+1]) )
				{	GenD0pT[k]++;} 
			}
			// D0 recosntruction Efficiency Eta
			double D0RangeEta[] = {0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.1};
			sizeArray = sizeof(D0RangeEta)/sizeof(D0RangeEta[0]); // Number of elements
			for( int k = 0; k < (sizeArray-1) ; k++)
			{	if( (MCD0Vec[0]->at(i) > D0RangeEta[k]) and (MCD0Vec[0]->at(i) < D0RangeEta[k+1]) )
				{	GenD0Eta[k]++;} 
			}

			t_D0MC.Fill();
			if (debug)cout << "debug D* 23 --------------------"<< endl;
			//--------------------------
			// D0 Matching
			if ( D0Vec[0]->size() > 0 ){
			double minDeltaR = 1000.;
			int id = -1; // Protection
			double deltaR = 0.;
			double deltaEta= 0.;
			double deltaPhi= 0.;
			for(unsigned int j=0; j < D0Vec[0]->size(); j++)
			{	
				//{"D0Kpimass", "D0Kpipt", "D0Kpieta", "D0Kpiphi", "D0lifetime", "D0Kpi_VtxProb", "D0KpiDispAngle", "D0KpisXY", "D0KpidXY", "D0Kpis3D", "D0Kpid3D", "D0KpikT"}			
				if( D0Vec[5]->at(j) < 0.01) continue;
				if( D0Vec[6]->at(j) < 0.99 ) continue;
				if( abs(D0Vec[0]->at(j)-1.86484) > 0.1) continue;
				if( D0Vec[1]->at(j) < 2. ) continue;
				if( D0Vec[9]->at(j) < SigCutD0) continue;

				//----------------------------------------
				//For D* Contamination			
				if (debug)cout << "For D0 Contamination"<< endl;
				if( FlagD0fromB->size() > 0)
				{	if( FlagD0fromB->at(i) == 1)
					{	D0Sec++;
						var_D0FromBmass = D0Vec[0]->at(j);
						var_D0FromBpt = D0Vec[1]->at(j);
						var_D0FromBeta = D0Vec[2]->at(j);
					}
					else
					{	D0Prompt++;
						var_D0FromPPmass = D0Vec[0]->at(j);
						var_D0FromPPpt = D0Vec[1]->at(j);
						var_D0FromPPeta = D0Vec[2]->at(j);
					}
				}
				t_D0Contamination.Fill();

				//-----------------------------------
				//DeltaR Function
				//{"MCpromptD0eta", "MCpromptD0phi", "MCpromptD0pt", "MCpromptD0energy", "MCpromptD0p", "MCpromptD0et", "MCpromptD0rapidity", "MCpromptD0mass", "MCpromptD0_DispAngle", "MCpromptD0lifetime" };
				deltaEta = MCD0Vec[0]->at(i) - D0Vec[2]->at(j);
				deltaPhi = MCD0Vec[1]->at(i) - D0Vec[3]->at(j);
				deltaR = sqrt( pow( deltaEta, 2) + pow( deltaPhi ,2));
				if( deltaR < minDeltaR) { id = j; minDeltaR = deltaR;} //Get the lower delta
				if (debug)cout << "debug D* 20 --------------------"<< endl;
			}

			if (id != -1){
				double recMass = -1; double GenMass = -1; // Protection
				double recPt = -1; double GenPt = -1;
				double recEta = -1; double GenEta = -1;
				double recTime = -1; double GenTime = -1;
				//Histograms for efficiency
				if( deltaR < 0.03){
					recMass = D0Vec[0]->at(id); GenMass = MCD0Vec[7]->at(i);
					recPt = D0Vec[1]->at(id); GenPt = MCD0Vec[2]->at(i);
					recEta = D0Vec[2]->at(id); GenEta = MCD0Vec[0]->at(i);
					recTime = ( (D0Vec[0]->at(id))*(D0Vec[8]->at(id))*pow(10,-2) ) / (c*(D0Vec[1]->at(id)) ); 
					GenTime = MCD0Vec[9]->at(i);

					hdeltaRD0Mass->Fill( recMass, GenMass );
					hdeltaRD0Pt->Fill( recPt, GenPt );
					hdeltaRD0Eta->Fill( abs(recEta), abs(GenEta) );
					hdeltaRD0Time->Fill( recTime, GenTime );

					// D* recosntruction Efficiency pT
					sizeArray = sizeof(D0RangePT)/sizeof(D0RangePT[0]); // Number of elements
					for( int k = 0; k < (sizeArray-1) ; k++)
					{	if( (D0Vec[1]->at(id) > D0RangePT[k] ) and ( D0Vec[1]->at(id) < D0RangePT[k+1]) )
						{	RecD0pT[k]++;} }
					// D* recosntruction Efficiency Eta
					sizeArray = sizeof(D0RangeEta)/sizeof(D0RangeEta[0]); // Number of elements
					for( int k = 0; k < (sizeArray-1) ; k++)
					{	if( (D0Vec[2]->at(id) > D0RangeEta[k] ) and (D0Vec[2]->at(id) < D0RangeEta[k+1]) )
						{	RecD0Eta[k]++;} }		
				}
		
			} // End ID
			} // End D* Matching


		} // End loop D0 MC
		} // End D0 MC protection
		
#endif

	if (debug)cout << "debug 24 --------------------"<< endl;	
	}//End loop tree entries for file f1


#if GENvaraibles == 0
#elif GENvaraibles == 1
	// Rec D* Efficiency
	TString DstarRangePTString[] = {"4", "5", "6", "7", "8", "12", "16", "24", "40", "100"};
	sizeArray = sizeof(DstarRangePTString)/sizeof(DstarRangePTString[0]); // Number of elements
	cout << "------Rec D* Efficiency--------------"<< endl;
	cout << "pT range / N Gen / N Rec / Efficiency  ----------"<< endl;
	for( int k = 0; k < (sizeArray-1) ; k++)
	{	cout << "$["+DstarRangePTString[k]+"-"+DstarRangePTString[k+1]+"]$ & " << GenDspT[k] << " & "<< RecDspT[k] << " & " << RecDspT[k]*1./GenDspT[k] << " end" << endl; }

	TString DstarRangeEtaString[] = {"0", "0.2", "0.4", "0.6", "0.8", "1.", "1.2", "1.4", "1.6", "1.8", "2.1"};
	sizeArray = sizeof(DstarRangeEtaString)/sizeof(DstarRangeEtaString[0]); // Number of elements
	cout << "eta range / N Gen / N Rec / Efficiency  ----------"<< endl;
	for( int k = 0; k < (sizeArray-1) ; k++)
	{	cout << "$["+DstarRangeEtaString[k]+"-"+DstarRangeEtaString[k+1]+"]$ & " << GenDsEta[k] << " & "<< RecDsEta[k] << " & " << RecDsEta[k]*1./GenDsEta[k] << " end" << endl; }

	// Rec D0 Efficiency
	TString D0RangePTString[] = {"4", "5", "6", "7", "8", "12", "16", "24", "40", "100"};
	sizeArray = sizeof(D0RangePTString)/sizeof(D0RangePTString[0]); // Number of elements
	cout << "------Rec D0 Efficiency--------------"<< endl;
	cout << "pT range / N Gen / N Rec / Efficiency  ----------"<< endl;
	for( int k = 0; k < (sizeArray-1) ; k++)
	{	cout << "$["+D0RangePTString[k]+"-"+D0RangePTString[k+1]+"]$ & " << GenD0pT[k] << " & "<< RecD0pT[k] << " & " << RecD0pT[k]*1./GenD0pT[k] << " end" << endl; }
	TString D0RangeEtaString[] = {"0", "0.2", "0.4", "0.6", "0.8", "1.", "1.2", "1.4", "1.6", "1.8", "2.1"};
	sizeArray = sizeof(D0RangeEtaString)/sizeof(D0RangeEtaString[0]); // Number of elements
	cout << "Eta range / N Gen / N Rec / Efficiency  ----------"<< endl;
	for( int k = 0; k < (sizeArray-1) ; k++)
	{	cout << "$["+D0RangeEtaString[k]+"-"+D0RangeEtaString[k+1]+"]$ & " << GenD0Eta[k] << " & "<< RecD0Eta[k] << " & " << RecD0Eta[k]*1./GenD0Eta[k] << " end" << endl; }
#endif


	TLatex* tex13;
	TCanvas* canvas = new TCanvas("PreliminarStudyDs","",900,600); canvas->Divide(3,3);

	//-----------------------------------------------------------------------------
	canvas = new TCanvas("PreliminarStudyDsMass","",900,600); canvas->Divide(3,3);	
	TString DsmassArray[] = { "= 0", "> 0.5", "> 1.0", ">1.5", "> 2.0", "> 2.5", "> 3.0"};
	sizeArray = sizeof(DsmassArray) / sizeof(DsmassArray[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ )
	{ 	canvas->cd(k+1);	hDstarMass[k]->Draw("e1p");
		tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig "+DsmassArray[k]+") }");
		if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig "+DsmassArray[k]+") }");}
		if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig "+DsmassArray[k]+") }");}
		tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();	 
	}
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMass_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMass_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 32 --------------------"<< endl;
	//-----------------------------------------------------------------------------
	canvas = new TCanvas("PreliminarStudyDsMassPt","",900,600); canvas->Divide(3,3);	
	TString DsmassPtArray[] = { "= 0", "> 0.5", "> 1.0", ">1.5", "> 2.0", "> 2.5", "> 3.0", "> 3.5", "> 4.0"};
	sizeArray = sizeof(DsmassPtArray) / sizeof(DsmassPtArray[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ )
	{ 	canvas->cd(k+1);	hDstarMassPt[k]->Draw("e1p");
		tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (pT "+DsmassPtArray[k]+") }");
		if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT "+DsmassPtArray[k]+") }");}
		if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT "+DsmassPtArray[k]+") }");}
		tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();
	}
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMassPt_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMassPt_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 33 --------------------"<< endl;
	//-----------------------------------------------------------------------------
	canvas = new TCanvas("PreliminarStudyD0Mass","",900,600); canvas->Divide(3,4);	
	TString D0KpimassArray[] = { "= 0", "> 0.5", "> 1.0", ">1.5", "> 2.0", "> 2.5", "> 3.0", "> 3.5", "> 4.0", "> 4.5", "> 5.0"};
	sizeArray = sizeof(hD0Mass) / sizeof(hD0Mass[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ )
	{ 	canvas->cd(k+1);	hD0Mass[k]->Draw("e1p");
		tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig "+D0KpimassArray[k]+") }");
		if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig "+D0KpimassArray[k]+") }");}
		if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT "+D0KpimassArray[k]+") }");}
		tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();	 }
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0Mass_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0Mass_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 34 --------------------"<< endl;
	//-----------------------------------------------------------------------------
	canvas = new TCanvas("PreliminarStudyD0MassPt","",900,600); canvas->Divide(3,3);
	sizeArray = sizeof(hD0MassPt) / sizeof(hD0MassPt[0]); // Number of elements 
	TString D0KpimassPtArray[] = { "= 0", "> 0.5", "> 1.0", ">1.5", "> 2.0", "> 2.5", "> 3.0", "> 3.5", "> 4.0" }; 
	for( int k=0; k < sizeArray; k++ )
	{ 	canvas->cd(k+1);	hD0MassPt[k]->Draw("e1p");
		tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (pT "+D0KpimassPtArray[k]+") }");
		if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT "+D0KpimassPtArray[k]+") }");}
		if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (pT "+D0KpimassPtArray[k]+") }");}
		tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();	 }
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0MassPt_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0MassPt_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));		
	if (debug)cout << "debug 35 --------------------"<< endl;
	//-----------------------------------------------------------------------------

	canvas = new TCanvas("StudyRecGenDstar","",900,600); canvas->Divide(2,2);
	TH1* nRecGenDstar[] = {hdeltaRDstarMass, hdeltaRDstarPt, hdeltaRDstarEta, hdeltaRDstarTime };
	sizeArray = sizeof(nRecGenDstar) / sizeof(nRecGenDstar[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ ){	
	canvas->cd(k+1);	nRecGenDstar[k]->Draw(); tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw(); }
	canvas->SaveAs(Form("%s%s/StudyRecGenDstar_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/StudyRecGenDstar_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 37 --------------------"<< endl;
	//**********************************************************
	canvas = new TCanvas("StudyRecGenD0","",900,600); canvas->Divide(2,2);
	TH1* nRecGenD0[] = {hdeltaRD0Mass, hdeltaRD0Pt, hdeltaRD0Eta, hdeltaRD0Time };
	sizeArray = sizeof(nRecGenD0) / sizeof(nRecGenD0[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ ){	
	canvas->cd(k+1);	nRecGenD0[k]->Draw(); tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.3");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw(); }
	canvas->SaveAs(Form("%s%s/StudyRecGenD0_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/StudyRecGenD0_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 38 --------------------"<< endl;

		cout << "test1" << endl;
	//if( TEfficiency::CheckConsistency( hRecEta_Rec1, hRecEta_MC ) )
	//{	//TCanvas* canvas = new TCanvas("canvas","",1200,600);
		//TGraphAsymmErrors *gr = new TGraphAsymmErrors( hRecEta_Rec1 , hRecEta_MC );
		//cout << "test2" << endl;
	//}


	f_analysis.cd();
	t_analysis.Write();  //Write in the root file
	t_D0analysis.Write();
	t_DsMC.Write();

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
TH2* makeTH2GenRec( TString name, TString nameAxis, Double_t BIN, Double_t NMIN , Double_t NMAX ){
	//Create a TH2D (one dimension) with alternatives bins
	TH2D *hh = new TH2D(name, name, BIN, NMIN, NMAX, BIN, NMIN, NMAX) ;
	hh->SetTitle(""+name+"; Rec "+nameAxis+"; Gen "+nameAxis+"" ); hh->SetMarkerStyle(1);
	hh->Sumw2(); hh->Sumw2(); return hh;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*void TGraphAsymmErrorsRec( const TH1 *h1, const TH1 *h2, const char* hTITLE, string path)
{	//Make a TGraphAsymmErrors using histograms as input

	if( TEfficiency::CheckConsistency(*h1,*h2) )
	{	TCanvas* canvas = new TCanvas("canvas","",1200,600);
		TGraphAsymmErrors *gr = new TGraphAsymmErrors( *(h1), *(h2));
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
}*/
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
