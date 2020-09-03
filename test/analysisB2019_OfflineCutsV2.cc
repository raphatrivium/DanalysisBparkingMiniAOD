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
#include "TGraphAsymmErrors.h"

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
TH1* makeTH1(const char* name, Double_t BIN, Double_t NMIN, Double_t NMAX, const char* TITLES );
TH1* makeTH1Rec(const char* name, const int nbins_low , double* nbins_high );
TH2* makeTH2GenRec( TString name, TString nameAxis, Double_t BIN, Double_t NMIN , Double_t NMAX );
void TGraphAsymmErrorsRec( TH1* h1, TH1* h2, const char* NAME, const char* hTITLE, int Nbins, string DATASETNAME, string path );
//--------------------------------------------------------------
// M A I N   F U N C T I O N 
//--------------------------------------------------------------
void analysisB2019_OfflineCutsV2()
{
	clock_t tStart = clock();

	//Datasets: BparkingDataset1 - MC_DStarToD0Pi_D0KPi - MC_MinBias - MC_BdToDStarX_ToD0Pi - MC_BuToDStarX_ToD0Pi
	string Dataset = "MC_DStarToD0Pi_D0KPi";
	//Put "" for select no trigger 
	string TriggerPath = ""; //HLT_Mu9_IP6
	if (Dataset != "MC_DStarToD0Pi_D0KPi" and Dataset != "MC_MinBias" and Dataset != "MC_BdToDStarX_ToD0Pi" and Dataset != "MC_BuToDStarX_ToD0Pi"  )
	{TriggerPath = "HLT_Mu9_IP6";} //HLT_Mu9_IP6

	//Save the file
	string path = "/eos/user/r/ragomesd/crab/haddteste/";
	string path2 = "/eos/user/r/ragomesd/analysisB2019/canvas/";

	double SigCutD0fromDs = 3. ; //3.
	double SigCutD0 = 5.; //5.

	cout << "SigCutD0fromDs: " << SigCutD0fromDs << endl;
	cout << "SigCutD0: " << SigCutD0 << endl;

	bool FlagDs = true;	bool FlagDsMC = true;
	bool FlagD0 = true;	bool FlagD0MC = true;
	bool debug = false;

	TChain* t1 = callTchain("./");
	Long64_t nentries = t1->GetEntries(); //Reading Number of tree entries of the file
	//nentries = 50000; //Test
	cout<< "Number of tree entries: "<< nentries << endl;
	Long64_t partpercent = nentries*0.05; //Percent done to show
	//-------Create root file-------------------------------------
	//TFile f_analysis("DMesonsHistograms.root","recreate"); //Creates root file
	TFile f_analysis(Form("%s%s.root", path.c_str(), Dataset.c_str() ),"recreate"); //Creates root file
	TTree t_analysis("t_analysis","t_analysis"); //Creates a Tree
	TTree t_DstarWrongCombination("t_DstarWrongCombination","t_DstarWrongCombination");
	TTree t_D0analysis("t_D0analysis","t_D0analysis");
	TTree t_D0WrongCombination("t_D0WrongCombination","t_D0WrongCombination");
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
	TString DstarString [] = {"Dsmass", "Dseta", "Dsphi", "Dspt", "DSDeltaR", "D0mass", "D0eta", "D0phi", "D0pt", "D0fromDSdXY", "D0fromDSsXY", "D0fromDSs3D", "D0fromDSd3D", "Anglephi", "D0_VtxProb", "Dslifetime", "DstarMinusD0mass", "TrkKpt", "TrkKeta", "TrkKphi", "Trkpipt", "Trkpieta", "Trkpiphi", "TrkSpt", "TrkSeta", "TrkSphi", };

	vector< vector<double>* > DstarVec;
	vector<double>* Dsmass =0; vector<double>* Dseta =0; vector<double>* Dsphi =0; vector<double>* Dspt =0; vector<double>* DSDeltaR =0; vector<double>* D0mass =0; vector<double>* D0eta =0; vector<double>* D0phi =0; vector<double>* D0pt =0; vector<double>* D0fromDSdXY =0; vector<double>* D0fromDSsXY =0; vector<double>* D0fromDSs3D =0; vector<double>* D0fromDSd3D =0; vector<double>* Anglephi =0; vector<double>* D0_VtxProb =0; vector<double>* Dslifetime =0; vector<double>* DstarMinusD0mass =0; vector<double>* TrkKpt =0; vector<double>* TrkKeta =0; vector<double>* TrkKphi =0; vector<double>* Trkpipt =0; vector<double>* Trkpieta =0; vector<double>* Trkpiphi =0; vector<double>* TrkSpt =0; vector<double>* TrkSeta =0; vector<double>* TrkSphi =0; 

	DstarVec.push_back(Dsmass); DstarVec.push_back(Dseta); DstarVec.push_back(Dsphi); DstarVec.push_back(Dspt); DstarVec.push_back(DSDeltaR); DstarVec.push_back(D0mass); DstarVec.push_back(D0eta); DstarVec.push_back(D0phi); DstarVec.push_back(D0pt); DstarVec.push_back(D0fromDSdXY); DstarVec.push_back(D0fromDSsXY); DstarVec.push_back(D0fromDSs3D); DstarVec.push_back(D0fromDSd3D); DstarVec.push_back(Anglephi); DstarVec.push_back(D0_VtxProb); DstarVec.push_back(Dslifetime); DstarVec.push_back(DstarMinusD0mass); DstarVec.push_back(TrkKpt); DstarVec.push_back(TrkKeta); DstarVec.push_back(TrkKphi); DstarVec.push_back(Trkpipt); DstarVec.push_back(Trkpieta); DstarVec.push_back(Trkpiphi); DstarVec.push_back(TrkSpt); DstarVec.push_back(TrkSeta); DstarVec.push_back(TrkSphi); 

	// Input Branches
	vector< TBranch* > DstarBranch; DstarBranch.clear();
	TBranch *b_Dsmass; TBranch *b_Dseta; TBranch *b_Dsphi; TBranch *b_Dspt; TBranch *b_DSDeltaR; TBranch *b_D0mass; TBranch *b_D0eta; TBranch *b_D0phi; TBranch *b_D0pt; TBranch *b_D0fromDSdXY; TBranch *b_D0fromDSsXY; TBranch *b_D0fromDSs3D; TBranch *b_D0fromDSd3D; TBranch *b_Anglephi; TBranch *b_D0_VtxProb; TBranch *b_Dslifetime; TBranch *b_DstarMinusD0mass; TBranch *b_TrkKpt; TBranch *b_TrkKeta; TBranch *b_TrkKphi; TBranch *b_Trkpipt; TBranch *b_Trkpieta; TBranch *b_Trkpiphi; TBranch *b_TrkSpt; TBranch *b_TrkSeta; TBranch *b_TrkSphi; 

	DstarBranch.push_back(b_Dsmass); DstarBranch.push_back(b_Dseta); DstarBranch.push_back(b_Dsphi); DstarBranch.push_back(b_Dspt); DstarBranch.push_back(b_DSDeltaR); DstarBranch.push_back(b_D0mass); DstarBranch.push_back(b_D0eta); DstarBranch.push_back(b_D0phi); DstarBranch.push_back(b_D0pt); DstarBranch.push_back(b_D0fromDSdXY); DstarBranch.push_back(b_D0fromDSsXY); DstarBranch.push_back(b_D0fromDSs3D); DstarBranch.push_back(b_D0fromDSd3D); DstarBranch.push_back(b_Anglephi); DstarBranch.push_back(b_D0_VtxProb); DstarBranch.push_back(b_Dslifetime); DstarBranch.push_back(b_DstarMinusD0mass); DstarBranch.push_back(b_TrkKpt); DstarBranch.push_back(b_TrkKeta); DstarBranch.push_back(b_TrkKphi); DstarBranch.push_back(b_Trkpipt); DstarBranch.push_back(b_Trkpieta); DstarBranch.push_back(b_Trkpiphi); DstarBranch.push_back(b_TrkSpt); DstarBranch.push_back(b_TrkSeta); DstarBranch.push_back(b_TrkSphi); 		
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
	double var_Dsmass =0.; double var_Dseta =0.; double var_Dsphi =0.; double var_Dspt =0.; double var_DSDeltaR =0.; double var_D0mass =0.; double var_D0eta =0.; double var_D0phi =0.; double var_D0pt =0.; double var_D0fromDSdXY =0.; double var_D0fromDSsXY =0.; double var_D0fromDSs3D =0.; double var_D0fromDSd3D =0.; double var_Anglephi =0.; double var_D0_VtxProb =0.; double var_Dslifetime =0.; double var_DstarMinusD0mass =0.; double var_TrkKpt =0.; double var_TrkKeta =0.; double var_TrkKphi =0.; double var_Trkpipt =0.; double var_Trkpieta =0.; double var_Trkpiphi =0.; double var_TrkSpt =0.; double var_TrkSeta =0.; double var_TrkSphi =0.;  
	double DstarVar[] = {var_Dsmass, var_Dseta, var_Dsphi, var_Dspt, var_DSDeltaR, var_D0mass, var_D0eta, var_D0phi, var_D0pt, var_D0fromDSdXY, var_D0fromDSsXY, var_D0fromDSs3D, var_D0fromDSd3D, var_Anglephi, var_D0_VtxProb, var_Dslifetime, var_DstarMinusD0mass, var_TrkKpt, var_TrkKeta, var_TrkKphi, var_Trkpipt, var_Trkpieta, var_Trkpiphi, var_TrkSpt, var_TrkSeta, var_TrkSphi, };

	sizeArray = sizeof(DstarVar)/sizeof(DstarVar[0]); //Number of elements
	vector<double> DstarVarVec; DstarVarVec.clear();
	DstarVarVec.insert(DstarVarVec.begin(), DstarVar, DstarVar+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < DstarVarVec.size(); k++)
	{ t_analysis.Branch( DstarString[k], &DstarVarVec[k], "var_"+DstarString[k]+"/D"); }

	double DstarCt = 0.;
	t_analysis.Branch("DstarCt",&DstarCt, "DstarCt/D");
#if GENvaraibles == 0	
#elif GENvaraibles == 1
	double PUWeight = 0.;
	TBranch *b_PUWeight; t1->SetBranchAddress("PUWeight",&PUWeight,&b_PUWeight);
	t_analysis.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif
	//----------------------------------------
	TString htitle = "Invariant Mass of the D* ; Mass [GeV] ; Events ";
	TH1 *DsmassHisto0 = makeTH1("DsmassHisto0", 100,1.93,2.1, htitle );
	TH1 *DsmassHisto1 = makeTH1("DsmassHisto1", 100,1.93,2.1, htitle );
	TH1 *DsmassHisto2 = makeTH1("DsmassHisto2", 100,1.93,2.1, htitle );
	TH1 *DsmassHisto3 = makeTH1("DsmassHisto3", 100,1.93,2.1, htitle );
	TH1 *DsmassHisto4 = makeTH1("DsmassHisto4", 100,1.93,2.1, htitle );
	TH1 *DsmassHisto5 = makeTH1("DsmassHisto5", 100,1.93,2.1, htitle );
	TH1 *DsmassHisto6 = makeTH1("DsmassHisto6", 100,1.93,2.1, htitle );
	TH1* hDstarMass[] = { DsmassHisto0, DsmassHisto1, DsmassHisto2, DsmassHisto3, DsmassHisto4, DsmassHisto5, DsmassHisto6};

	TH1 *DsmassPtHisto1 = makeTH1("DsmassPtHisto1", 100,1.93,2.1, htitle );
	TH1 *DsmassPtHisto2 = makeTH1("DsmassPtHisto2", 100,1.93,2.1, htitle );
	TH1 *DsmassPtHisto3 = makeTH1("DsmassPtHisto3", 100,1.93,2.1, htitle );
	TH1 *DsmassPtHisto4 = makeTH1("DsmassPtHisto4", 100,1.93,2.1, htitle );
	TH1 *DsmassPtHisto5 = makeTH1("DsmassPtHisto5", 100,1.93,2.1, htitle );
	TH1 *DsmassPtHisto6 = makeTH1("DsmassPtHisto6", 100,1.93,2.1, htitle );
	TH1 *DsmassPtHisto7 = makeTH1("DsmassPtHisto7", 100,1.93,2.1, htitle );
	TH1 *DsmassPtHisto8 = makeTH1("DsmassPtHisto8", 100,1.93,2.1, htitle );
	TH1 *DsmassPtHisto9 = makeTH1("DsmassPtHisto9", 100,1.93,2.1, htitle );
	TH1* hDstarMassPt[] = { DsmassPtHisto1, DsmassPtHisto2, DsmassPtHisto3, DsmassPtHisto4, DsmassPtHisto5, DsmassPtHisto6, DsmassPtHisto7, DsmassPtHisto8, DsmassPtHisto9};

	htitle = "LifeTime D0(from D*) ; t [s] ; Events ";
	TH1 *DsLifeTimeHisto0 = makeTH1("DsLifeTimeHisto0", 100, 0, 1.2*(pow(10,-12)), htitle );
	TH1 *DsLifeTimeHisto1 = makeTH1("DsLifeTimeHisto1", 100, 0, 1.2*(pow(10,-12)), htitle );
	TH1 *DsLifeTimeHisto2 = makeTH1("DsLifeTimeHisto2", 100, 0, 1.2*(pow(10,-12)), htitle );
	TH1 *DsLifeTimeHisto3 = makeTH1("DsLifeTimeHisto3", 100, 0, 1.2*(pow(10,-12)), htitle );
	TH1 *DsLifeTimeHisto4 = makeTH1("DsLifeTimeHisto4", 100, 0, 1.2*(pow(10,-12)), htitle );
	TH1 *DsLifeTimeHisto5 = makeTH1("DsLifeTimeHisto5", 100, 0, 1.2*(pow(10,-12)), htitle );
	TH1 *DsLifeTimeHisto6 = makeTH1("DsLifeTimeHisto6", 100, 0, 1.2*(pow(10,-12)), htitle );
	TH1* hDstarLifeTime[] = { DsLifeTimeHisto0, DsLifeTimeHisto1, DsLifeTimeHisto2, DsLifeTimeHisto3, DsLifeTimeHisto4, DsLifeTimeHisto5, DsLifeTimeHisto6};

	// Histograms D* 
	TH1 *hDsmass = makeTH1("hDsmass", 100, 1.93, 2.1, "Invariant Mass of the D* ; Mass [GeV/c2] ; Events " );
	TH1 *hDseta = makeTH1("hDseta", 100, -3., 3., "Pseudo-rapidity of the D0(From D*) ; #eta ; Events" );
	TH1 *hDsphi = makeTH1("hDsphi", 100, -4, 4, "Angle #phi of the D* ; #phi [rad] ; Events " );
	TH1 *hDspt = makeTH1("hDspt", 100, 2., 20., "Tranverse Momentum of the D* ; p_{T} [GeV] ; Events" );
	TH1 *hDSDeltaR = makeTH1("hDSDeltaR", 100, 0, 1., "#Delta R of the D* ; #Delta R ; Events" );
	TH1 *hhD0mass = makeTH1("hhD0mass", 100,1.76,1.96, "Invariant Mass of the D0(FromD*) ; Mass [GeV] ; Events " );
	TH1 *hD0eta = makeTH1("hD0eta", 100, -3., 3., "Invariant Mass of the D0(FromD*) ; #eta ; Events " );
	TH1 *hD0phi = makeTH1("hD0phi", 100, -4., 4., "Invariant Mass of the D0 ; #phi [rad] ; Events " );
	TH1 *hD0pt = makeTH1("hD0pt", 100, 0., 20., "Invariant Mass of the D0(FromD*) ; p_{T} [GeV] ; Events " );
	TH1 *hD0fromDSdXY = makeTH1("hD0fromDSdXY", 100, 0., 0.15, "Displacement of the D0(FromD*) ; Displacement XY [cm] ; Events " );
	TH1 *hD0fromDSsXY = makeTH1("hD0fromDSsXY", 100, 1., 6., "Significance (XY) of the D0(FromD*) ; Significance (XY) ; Events " );
	TH1 *hD0fromDSs3D = makeTH1("hD0fromDSs3D", 100, 1., 6., "Significance (3D) of the D0(FromD*) ; Significance (3D) ; Events " );
	TH1 *hD0fromDSd3D = makeTH1("hD0fromDSd3D", 100, 0., 0.25, "Displacement (3D) of the D0(FromD*) ; Displacement [cm] ; Events " );
	TH1 *hAnglephi = makeTH1("hAnglephi", 100, 0.989, 1.00, "Coseno Angle Phi of the D0(FromD*) ; cos#Phi ; Events " );
	TH1 *hD0_VtxProb = makeTH1("hD0_VtxProb", 100, 0, 1., "Vertex Probability of the D0(FromD*) ; Probability ; Events " );
	TH1 *hDslifetime = makeTH1("hDslifetime", 100, 0, 1.2*(pow(10,-12)), "Lifetime of the D0(FromD*) ; lifetime [s] ; Events " );
	TH1 *hDsMinusD0 = makeTH1("hDsMinusD0", 100,0.14,0.16, "#Delta m = m(K#pi#pi_{slow}) - m(K#pi) ; #Delta m [GeV/c2]; Events" );
	TH1 *hTrkKpt = makeTH1("hTrkKpt", 100, 0., 10., "Tranverse Momentum of Kaons ; p_{T} [GeV] ; Events" );
	TH1 *hTrkKeta = makeTH1("hTrkKeta", 100, -3., 3., "Pseudo-Rapidity of Kaons ; #eta ; Events" );
	TH1 *hTrkKphi = makeTH1("hTrkKphi", 100, -4, 4, "Angle #phi of Kaons ; #phi [rad] ; Events" );
	TH1 *hTrkpipt = makeTH1("hTrkpipt", 100, 0., 10., "Tranverse Momentum of Pions ; p_{T} [GeV] ; Events" );
	TH1 *hTrkpieta = makeTH1("hTrkpieta", 100, -3., 3., "Pseudo-Rapidity of Pions ; #eta ; Events" );
	TH1 *hTrkpiphi = makeTH1("hTrkpiphi", 100, -4, 4, "Angle #phi of Pions ; #phi [rad] ; Events" );
	TH1 *hTrkSpt = makeTH1("hTrkSpt", 100, 0., 4., "Tranverse Momentum of SlowPions ; p_{T} [GeV] ; Events" );
	TH1 *hTrkSeta = makeTH1("hTrkSeta", 100, -3., 3., "Pseudo-Rapidity of SlowPions ; #eta ; Events" );
	TH1 *hTrkSphi = makeTH1("hTrkSphi", 100, -4., 4., "Angle #phi of SlowPions ; #phi [rad] ; Events" );

	TH1* HistogramsDstar[] = { hDsmass, hDseta, hDsphi, hDspt, hDSDeltaR, hhD0mass, hD0eta , hD0phi, hD0pt, hD0fromDSdXY, hD0fromDSsXY, hD0fromDSs3D, hD0fromDSd3D, hAnglephi , hD0_VtxProb , hDslifetime , hDsMinusD0 , hTrkKpt , hTrkKeta , hTrkKphi , hTrkpipt , hTrkpieta , hTrkpiphi , hTrkSpt , hTrkSeta , hTrkSphi };


	//----------------------------------------
	//For D* Wrong combination
	//-----------------------------------------
	TString DstarWrongString [] = {"DsmassWrong", "DsetaWrong", "DsphiWrong", "DsptWrong", "D0massWrong", "D0etaWrong", "D0phiWrong", "D0ptWrong", "D0fromDSdXYWrong", "D0fromDSsXYWrong", "AnglephiWrong", "D0_VtxProbWrong", "DslifetimeWrong" };
	vector< vector<double>* > DstarVecWrong;
	vector<double>* DsmassWrong =0; vector<double>* DsetaWrong =0; vector<double>* DsphiWrong =0; vector<double>* DsptWrong =0; vector<double>* D0massWrong =0; vector<double>* D0etaWrong =0; vector<double>* D0phiWrong =0; vector<double>* D0ptWrong =0; vector<double>* D0fromDSdXYWrong =0; vector<double>* D0fromDSsXYWrong =0; vector<double>* AnglephiWrong =0; vector<double>* D0_VtxProbWrong =0; vector<double>* DslifetimeWrong =0; 
	DstarVecWrong.push_back(DsmassWrong); DstarVecWrong.push_back(DsetaWrong); DstarVecWrong.push_back(DsphiWrong); DstarVecWrong.push_back(DsptWrong); DstarVecWrong.push_back(D0massWrong); DstarVecWrong.push_back(D0etaWrong); DstarVecWrong.push_back(D0phiWrong); DstarVecWrong.push_back(D0ptWrong); DstarVecWrong.push_back(D0fromDSdXYWrong); DstarVecWrong.push_back(D0fromDSsXYWrong); DstarVecWrong.push_back(AnglephiWrong); DstarVecWrong.push_back(D0_VtxProbWrong); DstarVecWrong.push_back(DslifetimeWrong); 
	// Input Branches
	vector< TBranch* > DstarBranchWrong; DstarBranchWrong.clear();
	TBranch *b_DsmassWrong; TBranch *b_DsetaWrong; TBranch *b_DsphiWrong; TBranch *b_DsptWrong; TBranch *b_D0massWrong; TBranch *b_D0etaWrong; TBranch *b_D0phiWrong; TBranch *b_D0ptWrong; TBranch *b_D0fromDSdXYWrong; TBranch *b_D0fromDSsXYWrong; TBranch *b_AnglephiWrong; TBranch *b_D0_VtxProbWrong; TBranch *b_DslifetimeWrong; 
	DstarBranchWrong.push_back(b_DsmassWrong); DstarBranchWrong.push_back(b_DsetaWrong); DstarBranchWrong.push_back(b_DsphiWrong); DstarBranchWrong.push_back(b_DsptWrong); DstarBranchWrong.push_back(b_D0massWrong); DstarBranchWrong.push_back(b_D0etaWrong); DstarBranchWrong.push_back(b_D0phiWrong); DstarBranchWrong.push_back(b_D0ptWrong); DstarBranchWrong.push_back(b_D0fromDSdXYWrong); DstarBranchWrong.push_back(b_D0fromDSsXYWrong); DstarBranchWrong.push_back(b_AnglephiWrong); DstarBranchWrong.push_back(b_D0_VtxProbWrong); DstarBranchWrong.push_back(b_DslifetimeWrong);  
	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < DstarBranchWrong.size(); k++)
	{ t1->SetBranchAddress( DstarWrongString[k], &DstarVecWrong[k], &DstarBranchWrong[k]); }
	//-----------------------------------------
	//Variables 
	double var_DsmassWrong =0.; double var_DsetaWrong =0.; double var_DsphiWrong =0.; double var_DsptWrong =0.; double var_D0massWrong =0.; double var_D0etaWrong =0.; double var_D0phiWrong =0.; double var_D0ptWrong =0.; double var_D0fromDSdXYWrong =0.; double var_D0fromDSsXYWrong =0.; double var_AnglephiWrong =0.; double var_D0_VtxProbWrong =0.; double var_DslifetimeWrong =0.;  
	double DstarVarWrong[] = {var_DsmassWrong, var_DsetaWrong, var_DsphiWrong, var_DsptWrong, var_D0massWrong, var_D0etaWrong, var_D0phiWrong, var_D0ptWrong, var_D0fromDSdXYWrong, var_D0fromDSsXYWrong, var_AnglephiWrong, var_D0_VtxProbWrong, var_DslifetimeWrong };
	sizeArray = sizeof(DstarVarWrong)/sizeof(DstarVarWrong[0]); //Number of elements
	vector<double> DstarVarVecWrong; DstarVarVecWrong.clear();
	DstarVarVecWrong.insert(DstarVarVecWrong.begin(), DstarVarWrong, DstarVarWrong+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < DstarVarVecWrong.size(); k++)
	{ t_DstarWrongCombination.Branch( DstarWrongString[k], &DstarVarVecWrong[k], "var_"+DstarWrongString[k]+"/D"); }

	double DstarCtWrong =0.;
	t_DstarWrongCombination.Branch("DstarCtWrong",&DstarCtWrong, "DstarCtWrong/D");

#if GENvaraibles == 0	
#elif GENvaraibles == 1
	t_DstarWrongCombination.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif


	// Histograms D* Wrong 
	TH1 *hDsmassWrong = makeTH1("hDsmassWrong", 100, 1.93, 2.1, "Invariant Mass of the D* ; Mass [GeV/c2] ; Events " );
	TH1 *hDsetaWrong = makeTH1("hDsetaWrong", 100, -3., 3., "Pseudo-rapidity of the D0(From D*) ; #eta ; Events" );
	TH1 *hDsphiWrong = makeTH1("hDsphiWrong", 100, -4, 4, "Angle #phi of the D* ; #phi [rad] ; Events " );
	TH1 *hDsptWrong = makeTH1("hDsptWrong", 100, 2., 20., "Tranverse Momentum of the D* ; p_{T} [GeV] ; Events" );
	TH1 *hD0massWrong = makeTH1("hD0massWrong", 100,1.76,1.96, "Invariant Mass of the D0(FromD*) ; Mass [GeV] ; Events " );
	TH1 *hD0etaWrong = makeTH1("hD0etaWrong", 100, -3., 3., "Invariant Mass of the D0(FromD*) ; #eta ; Events " );
	TH1 *hD0phiWrong = makeTH1("hD0phiWrong", 100, -4., 4., "Invariant Mass of the D0 ; #phi [rad] ; Events " );
	TH1 *hD0ptWrong = makeTH1("hD0ptWrong", 100, 0., 20., "Invariant Mass of the D0(FromD*) ; p_{T} [GeV] ; Events " );
	TH1 *hD0fromDSdXYWrong = makeTH1("hD0fromDSdXYWrong", 100, 0., 0.15, "Displacement of the D0(FromD*) ; Displacement XY [cm] ; Events " );
	TH1 *hD0fromDSsXYWrong = makeTH1("hD0fromDSsXYWrong", 100, 1., 6., "Significance (XY) of the D0(FromD*) ; Significance (XY) ; Events " );
	TH1 *hAnglephiWrong = makeTH1("hAnglephiWrong", 100, 0.989, 1.00, "Coseno Angle Phi of the D0(FromD*) ; cos#Phi ; Events " );
	TH1 *hD0_VtxProbWrong = makeTH1("hD0_VtxProbWrong", 100, 0, 1., "Vertex Probability of the D0(FromD*) ; Probability ; Events " );
	TH1 *hDslifetimeWrong = makeTH1("hDslifetimeWrong", 100, 0, 1.2*(pow(10,-12)), "Lifetime of the D0(FromD*) ; lifetime [s] ; Events " );

	TH1* HistogramsDstarWrong[] = { hDsmassWrong, hDsetaWrong, hDsphiWrong, hDsptWrong, hD0massWrong, hD0etaWrong, hD0phiWrong , hD0ptWrong, hD0fromDSdXYWrong, hD0fromDSsXYWrong, hAnglephiWrong, hD0_VtxProbWrong, hDslifetimeWrong};
	

	//----------------------------------------
	//For D0
	//-----------------------------------------
	TString D0String [] = {"D0Kpimass", "D0Kpipt", "D0Kpieta", "D0Kpiphi", "D0lifetime", "D0Kpi_VtxProb", "D0KpiDispAngle", "D0KpisXY", "D0KpidXY", "D0Kpis3D", "D0Kpid3D", "D0KpikT", "TrkD0Kpt", "TrkD0Keta", "TrkD0kphi", "TrkD0pipt", "TrkD0pieta", "TrkD0piphi", };

	vector< vector<double>* > D0Vec;
	vector<double>* D0Kpimass =0; vector<double>* D0Kpipt =0; vector<double>* D0Kpieta =0; vector<double>* D0Kpiphi =0; vector<double>* D0lifetime =0; vector<double>* D0Kpi_VtxProb =0; vector<double>* D0KpiDispAngle =0; vector<double>* D0KpisXY =0; vector<double>* D0KpidXY =0; vector<double>* D0Kpis3D =0; vector<double>* D0Kpid3D =0; vector<double>* D0KpikT =0; vector<double>* TrkD0Kpt =0; vector<double>* TrkD0Keta =0; vector<double>* TrkD0kphi =0; vector<double>* TrkD0pipt =0; vector<double>* TrkD0pieta =0; vector<double>* TrkD0piphi =0; 

	D0Vec.push_back(D0Kpimass); D0Vec.push_back(D0Kpipt); D0Vec.push_back(D0Kpieta); D0Vec.push_back(D0Kpiphi); D0Vec.push_back(D0lifetime); D0Vec.push_back(D0Kpi_VtxProb); D0Vec.push_back(D0KpiDispAngle); D0Vec.push_back(D0KpisXY); D0Vec.push_back(D0KpidXY); D0Vec.push_back(D0Kpis3D); D0Vec.push_back(D0Kpid3D); D0Vec.push_back(D0KpikT); D0Vec.push_back(TrkD0Kpt); D0Vec.push_back(TrkD0Keta); D0Vec.push_back(TrkD0kphi); D0Vec.push_back(TrkD0pipt); D0Vec.push_back(TrkD0pieta); D0Vec.push_back(TrkD0piphi); 

	// Input Branches
	vector< TBranch* > D0Branch; D0Branch.clear();
	TBranch *b_D0Kpimass; TBranch *b_D0Kpipt; TBranch *b_D0Kpieta; TBranch *b_D0Kpiphi; TBranch *b_D0lifetime; TBranch *b_D0Kpi_VtxProb; TBranch *b_D0KpiDispAngle; TBranch *b_D0KpisXY; TBranch *b_D0KpidXY; TBranch *b_D0Kpis3D; TBranch *b_D0Kpid3D; TBranch *b_D0KpikT; TBranch *b_TrkD0Kpt; TBranch *b_TrkD0Keta; TBranch *b_TrkD0kphi; TBranch *b_TrkD0pipt; TBranch *b_TrkD0pieta; TBranch *b_TrkD0piphi; 

	D0Branch.push_back(b_D0Kpimass); D0Branch.push_back(b_D0Kpipt); D0Branch.push_back(b_D0Kpieta); D0Branch.push_back(b_D0Kpiphi); D0Branch.push_back(b_D0lifetime); D0Branch.push_back(b_D0Kpi_VtxProb); D0Branch.push_back(b_D0KpiDispAngle); D0Branch.push_back(b_D0KpisXY); D0Branch.push_back(b_D0KpidXY); D0Branch.push_back(b_D0Kpis3D); D0Branch.push_back(b_D0Kpid3D); D0Branch.push_back(b_D0KpikT); D0Branch.push_back(b_TrkD0Kpt); D0Branch.push_back(b_TrkD0Keta); D0Branch.push_back(b_TrkD0kphi); D0Branch.push_back(b_TrkD0pipt); D0Branch.push_back(b_TrkD0pieta); D0Branch.push_back(b_TrkD0piphi);
	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < D0Branch.size(); k++)
	{ t1->SetBranchAddress( D0String[k], &D0Vec[k], &D0Branch[k]); }
	//-----------------------------------------
	//Variables 
	double var_D0Kpimass =0.; double var_D0Kpipt =0.; double var_D0Kpieta =0.; double var_D0Kpiphi =0.; double var_D0lifetime =0.; double var_D0Kpi_VtxProb =0.; double var_D0KpiDispAngle =0.; double var_D0KpisXY =0.; double var_D0KpidXY =0.; double var_D0Kpis3D =0.; double var_D0Kpid3D =0.; double var_D0KpikT =0.; double var_TrkD0Kpt =0.; double var_TrkD0Keta =0.; double var_TrkD0kphi =0.; double var_TrkD0pipt =0.; double var_TrkD0pieta =0.; double var_TrkD0piphi =0.;  
	double D0Var[] = {var_D0Kpimass, var_D0Kpipt, var_D0Kpieta, var_D0Kpiphi, var_D0lifetime, var_D0Kpi_VtxProb, var_D0KpiDispAngle, var_D0KpisXY, var_D0KpidXY, var_D0Kpis3D, var_D0Kpid3D, var_D0KpikT, var_TrkD0Kpt, var_TrkD0Keta, var_TrkD0kphi, var_TrkD0pipt, var_TrkD0pieta, var_TrkD0piphi, };

	sizeArray = sizeof(D0Var)/sizeof(D0Var[0]); //Number of elements
	vector<double> D0VarVec; D0VarVec.clear();
	D0VarVec.insert(D0VarVec.begin(), D0Var, D0Var+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < D0VarVec.size(); k++)
	{ t_D0analysis.Branch( D0String[k], &D0VarVec[k], "var_"+D0String[k]+"/D"); }

	double D0Ct =0.;
	t_D0analysis.Branch("D0Ct",&D0Ct, "D0Ct/D");

#if GENvaraibles == 0	
#elif GENvaraibles == 1
	t_D0analysis.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif

	// Histograms D0 
	TH1 *hD0Kpimass = makeTH1("hD0Kpimass", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events " );
	TH1 *hD0Kpipt = makeTH1("hD0Kpipt", 100, 0., 20., "Invariant Mass of the D0 ; p_{T} [GeV] ; Events " );
	TH1 *hD0Kpieta = makeTH1("hD0Kpieta", 100, -3., 3., "Invariant Mass of the D0 ; #eta [GeV] ; Events " );
	TH1 *hD0Kpiphi = makeTH1("hD0Kpiphi", 100, -4., 4., "Invariant Mass of the D0 ; #phi [rad] ; Events " );
	TH1 *hD0lifetime = makeTH1("hD0lifetime", 100, 0, 2*(pow(10,-12)), "Lifetime of the D0 ; lifetime [s] ; Events " );
	TH1 *hD0Kpi_VtxProb = makeTH1("hD0Kpi_VtxProb", 100, 0, 1., "Vertex Probability of the D0 ; Probability ; Events " );
	TH1 *hD0KpiDispAngle = makeTH1("hD0KpiDispAngle", 100, 0.989, 1.00, "Coseno Angle Phi of the D0 ; cos#Phi ; Events " );
	TH1 *hD0KpisXY = makeTH1("hD0KpisXY", 100, 0., 20., "Significance (XY) of the D0 ; Significance ; Events " );
	TH1 *hD0KpidXY = makeTH1("hD0KpidXY", 100, 0., 0.15, "Displacement of the D0 ; Displacement XY [cm] ; Events " );
	TH1 *hD0Kpis3D = makeTH1("hD0Kpis3D", 100, 0., 20., "Significance (3D) of the D0 ; Significance (3D) ; Events " );
	TH1 *hD0Kpid3D = makeTH1("hD0Kpid3D", 100, 0., 0.5, "Displacement (3D) of the D0 ; Displacement [cm] ; Events " );
	TH1 *hD0KpikT = makeTH1("hD0KpikT", 100, 0.54, 1.0, "Vetorial product of Kaon and D0 ; K #times D0 ; Events" );
	TH1 *hTrkD0Kpt = makeTH1("hTrkD0Kpt", 100, 0., 10., "Tranverse Momentum of Kaons ; p_{T} [GeV] ; Events" );
	TH1 *hTrkD0Keta = makeTH1("hTrkD0Keta", 100, -3., 3., "Pseudo-Rapidity of Kaons ; #eta ; Events" );
	TH1 *hTrkD0kphi = makeTH1("hTrkD0kphi", 100, -4, 4, "Angle #phi of Kaons ; #phi [rad] ; Events" );
	TH1 *hTrkD0pipt = makeTH1("hTrkD0pipt", 100, 0., 10., "Tranverse Momentum of Pions ; p_{T} [GeV] ; Events" );
	TH1 *hTrkD0pieta = makeTH1("hTrkD0pieta", 100, -3., 3., "Pseudo-Rapidity of Pions ; #eta ; Events" );
	TH1 *hTrkD0piphi = makeTH1("hTrkD0piphi", 100, -4, 4, "Angle #phi of Pions ; #phi [rad] ; Events" );

	TH1* HistogramsD0[] = {hD0Kpimass, hD0Kpipt, hD0Kpieta, hD0Kpiphi, hD0lifetime, hD0Kpi_VtxProb, hD0KpiDispAngle , hD0KpisXY, hD0KpidXY, hD0Kpis3D, hD0Kpid3D, hD0KpikT, hTrkD0Kpt, hTrkD0Keta, hTrkD0kphi, hTrkD0pipt, hTrkD0pieta, hTrkD0piphi };

	//Study of Mass with different D0 significance cuts
	htitle = "Invariant Mass of the D0 ; Mass [GeV] ; Events ";
	TH1 *D0KpimassHisto0 = makeTH1("D0KpimassHisto0", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto1 = makeTH1("D0KpimassHisto1", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto2 = makeTH1("D0KpimassHisto2", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto3 = makeTH1("D0KpimassHisto3", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto4 = makeTH1("D0KpimassHisto4", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto5 = makeTH1("D0KpimassHisto5", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto6 = makeTH1("D0KpimassHisto6", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto7 = makeTH1("D0KpimassHisto7", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto8 = makeTH1("D0KpimassHisto8", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto9 = makeTH1("D0KpimassHisto9", 100,1.76,1.96, htitle );
	TH1 *D0KpimassHisto10 = makeTH1("D0KpimassHisto10", 100,1.76,1.96, htitle );
	TH1* hD0Mass[] = { D0KpimassHisto0, D0KpimassHisto1, D0KpimassHisto2, D0KpimassHisto3, D0KpimassHisto4, D0KpimassHisto5, D0KpimassHisto6, D0KpimassHisto7, D0KpimassHisto8, D0KpimassHisto9, D0KpimassHisto10};

	TH1 *D0KpimassPtHisto1 = makeTH1("D0KpimassPtHisto1", 100,1.76,1.96, htitle );
	TH1 *D0KpimassPtHisto2 = makeTH1("D0KpimassPtHisto2", 100,1.76,1.96, htitle );
	TH1 *D0KpimassPtHisto3 = makeTH1("D0KpimassPtHisto3", 100,1.76,1.96, htitle );
	TH1 *D0KpimassPtHisto4 = makeTH1("D0KpimassPtHisto4", 100,1.76,1.96, htitle );
	TH1 *D0KpimassPtHisto5 = makeTH1("D0KpimassPtHisto5", 100,1.76,1.96, htitle );
	TH1 *D0KpimassPtHisto6 = makeTH1("D0KpimassPtHisto6", 100,1.76,1.96, htitle );
	TH1 *D0KpimassPtHisto7 = makeTH1("D0KpimassPtHisto7", 100,1.76,1.96, htitle );
	TH1 *D0KpimassPtHisto8 = makeTH1("D0KpimassPtHisto8", 100,1.76,1.96, htitle );
	TH1 *D0KpimassPtHisto9 = makeTH1("D0KpimassPtHisto9", 100,1.76,1.96, htitle );
	TH1* hD0MassPt[] = { D0KpimassPtHisto1, D0KpimassPtHisto2, D0KpimassPtHisto3, D0KpimassPtHisto4, D0KpimassPtHisto5, D0KpimassPtHisto6, D0KpimassPtHisto7, D0KpimassPtHisto8, D0KpimassPtHisto9 };

	//----------------------------------------
	//For D0 Wrong Combination
	//-----------------------------------------
	TString D0StringWrong [] = {"D0KpimassWrong", "D0KpiptWrong", "D0KpietaWrong", "D0KpiphiWrong", "D0lifetimeWrong", "D0Kpi_VtxProbWrong", "D0KpiDispAngleWrong", "D0KpisXYWrong", "D0KpidXYWrong", "D0Kpis3DWrong", "D0Kpid3DWrong", "D0KpikTWrong" };

	vector< vector<double>* > D0VecWrong;
	vector<double>* D0KpimassWrong =0; vector<double>* D0KpiptWrong =0; vector<double>* D0KpietaWrong =0; vector<double>* D0KpiphiWrong =0; vector<double>* D0lifetimeWrong =0; vector<double>* D0Kpi_VtxProbWrong =0; vector<double>* D0KpiDispAngleWrong =0; vector<double>* D0KpisXYWrong =0; vector<double>* D0KpidXYWrong =0; vector<double>* D0Kpis3DWrong =0; vector<double>* D0Kpid3DWrong =0; vector<double>* D0KpikTWrong =0; 

	D0VecWrong.push_back(D0KpimassWrong); D0VecWrong.push_back(D0KpiptWrong); D0VecWrong.push_back(D0KpietaWrong); D0VecWrong.push_back(D0KpiphiWrong); D0VecWrong.push_back(D0lifetimeWrong); D0VecWrong.push_back(D0Kpi_VtxProbWrong); D0VecWrong.push_back(D0KpiDispAngleWrong); D0VecWrong.push_back(D0KpisXYWrong); D0VecWrong.push_back(D0KpidXYWrong); D0VecWrong.push_back(D0Kpis3DWrong); D0VecWrong.push_back(D0Kpid3DWrong); D0VecWrong.push_back(D0KpikTWrong); 

	// Input Branches
	vector< TBranch* > D0BranchWrong; D0BranchWrong.clear();
	TBranch *b_D0KpimassWrong; TBranch *b_D0KpiptWrong; TBranch *b_D0KpietaWrong; TBranch *b_D0KpiphiWrong; TBranch *b_D0lifetimeWrong; TBranch *b_D0Kpi_VtxProbWrong; TBranch *b_D0KpiDispAngleWrong; TBranch *b_D0KpisXYWrong; TBranch *b_D0KpidXYWrong; TBranch *b_D0Kpis3DWrong; TBranch *b_D0Kpid3DWrong; TBranch *b_D0KpikTWrong; 

	D0BranchWrong.push_back(b_D0KpimassWrong); D0BranchWrong.push_back(b_D0KpiptWrong); D0BranchWrong.push_back(b_D0KpietaWrong); D0BranchWrong.push_back(b_D0KpiphiWrong); D0BranchWrong.push_back(b_D0lifetimeWrong); D0BranchWrong.push_back(b_D0Kpi_VtxProbWrong); D0BranchWrong.push_back(b_D0KpiDispAngleWrong); D0BranchWrong.push_back(b_D0KpisXYWrong); D0BranchWrong.push_back(b_D0KpidXYWrong); D0BranchWrong.push_back(b_D0Kpis3DWrong); D0BranchWrong.push_back(b_D0Kpid3DWrong); D0BranchWrong.push_back(b_D0KpikTWrong); 

	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < D0BranchWrong.size(); k++)
	{ t1->SetBranchAddress( D0StringWrong[k], &D0VecWrong[k], &D0BranchWrong[k]); }
	//-----------------------------------------
	//Variables 
	double var_D0KpimassWrong =0.; double var_D0KpiptWrong =0.; double var_D0KpietaWrong =0.; double var_D0KpiphiWrong =0.; double var_D0lifetimeWrong =0.; double var_D0Kpi_VtxProbWrong =0.; double var_D0KpiDispAngleWrong =0.; double var_D0KpisXYWrong =0.; double var_D0KpidXYWrong =0.; double var_D0Kpis3DWrong =0.; double var_D0Kpid3DWrong =0.; double var_D0KpikTWrong =0.;  
	double D0VarWrong[] = {var_D0KpimassWrong, var_D0KpiptWrong, var_D0KpietaWrong, var_D0KpiphiWrong, var_D0lifetimeWrong, var_D0Kpi_VtxProbWrong, var_D0KpiDispAngleWrong, var_D0KpisXYWrong, var_D0KpidXYWrong, var_D0Kpis3DWrong, var_D0Kpid3DWrong, var_D0KpikTWrong };

	sizeArray = sizeof(D0VarWrong)/sizeof(D0VarWrong[0]); //Number of elements
	vector<double> D0VarVecWrong; D0VarVecWrong.clear();
	D0VarVecWrong.insert(D0VarVecWrong.begin(), D0VarWrong, D0VarWrong+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < D0VarVecWrong.size(); k++)
	{ t_D0WrongCombination.Branch( D0StringWrong[k], &D0VarVecWrong[k], "var_"+D0StringWrong[k]+"/D"); }

	double D0CtWrong =0.;
	t_D0WrongCombination.Branch("D0CtWrong",&D0CtWrong, "D0CtWrong/D");

#if GENvaraibles == 0	
#elif GENvaraibles == 1
	t_D0WrongCombination.Branch("PUWeight",&PUWeight, "PUWeight/D");
#endif

	// Histograms D0 
	TH1 *hD0KpimassWrong = makeTH1("hD0KpimassWrong", 100,1.76,1.96, "Invariant Mass of the D0 ; Mass [GeV] ; Events " );
	TH1 *hD0KpiptWrong = makeTH1("hD0KpiptWrong", 100, 0., 20., "Invariant Mass of the D0 ; p_{T} [GeV] ; Events " );
	TH1 *hD0KpietaWrong = makeTH1("hD0KpietaWrong", 100, -3., 3., "Invariant Mass of the D0 ; #eta [GeV] ; Events " );
	TH1 *hD0KpiphiWrong = makeTH1("hD0KpiphiWrong", 100, -4., 4., "Invariant Mass of the D0 ; #phi [rad] ; Events " );
	TH1 *hD0lifetimeWrong = makeTH1("hD0lifetimeWrong", 100, 0, 2*(pow(10,-12)), "Lifetime of the D0 ; lifetime [s] ; Events " );
	TH1 *hD0Kpi_VtxProbWrong = makeTH1("hD0Kpi_VtxProbWrong", 100, 0, 1., "Vertex Probability of the D0 ; Probability ; Events " );
	TH1 *hD0KpiDispAngleWrong = makeTH1("hD0KpiDispAngleWrong", 100, 0.989, 1.00, "Coseno Angle Phi of the D0 ; cos#Phi ; Events " );
	TH1 *hD0KpisXYWrong = makeTH1("hD0KpisXYWrong", 100, 0., 20., "Significance (XY) of the D0 ; Significance ; Events " );
	TH1 *hD0KpidXYWrong = makeTH1("hD0KpidXYWrong", 100, 0., 0.15, "Displacement of the D0 ; Displacement XY [cm] ; Events " );
	TH1 *hD0Kpis3DWrong = makeTH1("hD0Kpis3DWrong", 100, 0., 20., "Significance (3D) of the D0 ; Significance (3D) ; Events " );
	TH1 *hD0Kpid3DWrong = makeTH1("hD0Kpid3DWrong", 100, 0., 0.5, "Displacement (3D) of the D0 ; Displacement [cm] ; Events " );
	TH1 *hD0KpikTWrong = makeTH1("hD0KpikTWrong", 100, 0.54, 1.0, "Vetorial product of Kaon and D0 ; K #times D0 ; Events" );

	TH1* HistogramsD0Wrong[] = {hD0KpimassWrong, hD0KpiptWrong, hD0KpietaWrong, hD0KpiphiWrong, hD0lifetimeWrong, hD0Kpi_VtxProbWrong, hD0KpiDispAngleWrong , hD0KpisXYWrong, hD0KpidXYWrong, hD0Kpis3DWrong, hD0Kpid3DWrong, hD0KpikTWrong };


#if GENvaraibles == 0	
#elif GENvaraibles == 1
	//----------------------------------------
	//For D* MC
	//-----------------------------------------
		TString MCDstarString [] = {"MCDseta", "MCDsphi", "MCDspt", "MCDsenergy", "MCDsp", "MCDset", "MCDsmass", "MCD0eta", "MCD0phi", "MCD0pt", "MCD0energy", "MCD0p", "MCD0et", "MCD0rapidity", "MCD0mass", "MCD0lifetime", "MCD0dispXY", "MCDsSeta", "MCDsSpt", "MCDsKeta", "MCDsKpt", "MCDsPieta", "MCDsPipt" };

	vector< vector<double>* > MCDstarVec;
	vector<double>* MCDseta =0; vector<double>* MCDsphi =0; vector<double>* MCDspt =0; vector<double>* MCDsenergy =0; vector<double>* MCDsp =0; vector<double>* MCDset =0; vector<double>* MCDsmass =0; vector<double>* MCD0eta =0; vector<double>* MCD0phi =0; vector<double>* MCD0pt =0; vector<double>* MCD0energy =0; vector<double>* MCD0p =0; vector<double>* MCD0et =0; vector<double>* MCD0rapidity =0; vector<double>* MCD0mass =0; vector<double>* MCD0lifetime =0; vector<double>* MCD0dispXY =0; vector<double>* MCDsSeta =0; vector<double>* MCDsSpt =0; vector<double>* MCDsKeta =0; vector<double>* MCDsKpt =0; vector<double>* MCDsPieta =0; vector<double>* MCDsPipt =0; 

	MCDstarVec.push_back(MCDseta); MCDstarVec.push_back(MCDsphi); MCDstarVec.push_back(MCDspt); MCDstarVec.push_back(MCDsenergy); MCDstarVec.push_back(MCDsp); MCDstarVec.push_back(MCDset); MCDstarVec.push_back(MCDsmass); MCDstarVec.push_back(MCD0eta); MCDstarVec.push_back(MCD0phi); MCDstarVec.push_back(MCD0pt); MCDstarVec.push_back(MCD0energy); MCDstarVec.push_back(MCD0p); MCDstarVec.push_back(MCD0et); MCDstarVec.push_back(MCD0rapidity); MCDstarVec.push_back(MCD0mass); MCDstarVec.push_back(MCD0lifetime); MCDstarVec.push_back(MCD0dispXY); MCDstarVec.push_back(MCDsSeta); MCDstarVec.push_back(MCDsSpt); MCDstarVec.push_back(MCDsKeta); MCDstarVec.push_back(MCDsKpt); MCDstarVec.push_back(MCDsPieta); MCDstarVec.push_back(MCDsPipt); 

	// Input Branches
	vector< TBranch* > MCDstarBranch; MCDstarBranch.clear();
	TBranch *b_MCDseta; TBranch *b_MCDsphi; TBranch *b_MCDspt; TBranch *b_MCDsenergy; TBranch *b_MCDsp; TBranch *b_MCDset; TBranch *b_MCDsmass; TBranch *b_MCD0eta; TBranch *b_MCD0phi; TBranch *b_MCD0pt; TBranch *b_MCD0energy; TBranch *b_MCD0p; TBranch *b_MCD0et; TBranch *b_MCD0rapidity; TBranch *b_MCD0mass; TBranch *b_MCD0lifetime; TBranch *b_MCD0dispXY; TBranch *b_MCDsSeta; TBranch *b_MCDsSpt; TBranch *b_MCDsKeta; TBranch *b_MCDsKpt; TBranch *b_MCDsPieta; TBranch *b_MCDsPipt; 

	MCDstarBranch.push_back(b_MCDseta); MCDstarBranch.push_back(b_MCDsphi); MCDstarBranch.push_back(b_MCDspt); MCDstarBranch.push_back(b_MCDsenergy); MCDstarBranch.push_back(b_MCDsp); MCDstarBranch.push_back(b_MCDset); MCDstarBranch.push_back(b_MCDsmass); MCDstarBranch.push_back(b_MCD0eta); MCDstarBranch.push_back(b_MCD0phi); MCDstarBranch.push_back(b_MCD0pt); MCDstarBranch.push_back(b_MCD0energy); MCDstarBranch.push_back(b_MCD0p); MCDstarBranch.push_back(b_MCD0et); MCDstarBranch.push_back(b_MCD0rapidity); MCDstarBranch.push_back(b_MCD0mass); MCDstarBranch.push_back(b_MCD0lifetime); MCDstarBranch.push_back(b_MCD0dispXY); MCDstarBranch.push_back(b_MCDsSeta); MCDstarBranch.push_back(b_MCDsSpt); MCDstarBranch.push_back(b_MCDsKeta); MCDstarBranch.push_back(b_MCDsKpt); MCDstarBranch.push_back(b_MCDsPieta); MCDstarBranch.push_back(b_MCDsPipt); 
	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < MCDstarBranch.size(); k++)
	{ t1->SetBranchAddress( MCDstarString[k], &MCDstarVec[k], &MCDstarBranch[k]); }

	//-----------------------------------------
	//Variables 
	double var_MCDseta =0.; double var_MCDsphi =0.; double var_MCDspt =0.; double var_MCDsenergy =0.; double var_MCDsp =0.; double var_MCDset =0.; double var_MCDsmass =0.; double var_MCD0eta =0.; double var_MCD0phi =0.; double var_MCD0pt =0.; double var_MCD0energy =0.; double var_MCD0p =0.; double var_MCD0et =0.; double var_MCD0rapidity =0.; double var_MCD0mass =0.; double var_MCD0lifetime =0.; double var_MCD0dispXY =0.; double var_MCDsSeta =0.; double var_MCDsSpt =0.; double var_MCDsKeta =0.; double var_MCDsKpt =0.; double var_MCDsPieta =0.; double var_MCDsPipt =0.;  
	double MCDstarVar[] = {var_MCDseta, var_MCDsphi, var_MCDspt, var_MCDsenergy, var_MCDsp, var_MCDset, var_MCDsmass, var_MCD0eta, var_MCD0phi, var_MCD0pt, var_MCD0energy, var_MCD0p, var_MCD0et, var_MCD0rapidity, var_MCD0mass, var_MCD0lifetime, var_MCD0dispXY, var_MCDsSeta, var_MCDsSpt, var_MCDsKeta, var_MCDsKpt, var_MCDsPieta, var_MCDsPipt };

	sizeArray = sizeof(MCDstarVar)/sizeof(MCDstarVar[0]); //Number of elements
	vector<double> MCDstarVarVec; MCDstarVarVec.clear();
	MCDstarVarVec.insert(MCDstarVarVec.begin(), MCDstarVar, MCDstarVar+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < MCDstarVarVec.size(); k++)
	{ t_DsMC.Branch( MCDstarString[k], &MCDstarVarVec[k], "var_"+MCDstarString[k]+"/D"); }
	t_DsMC.Branch("PUWeight",&PUWeight, "PUWeight/D");

	double MCDstarCt =0.;
	t_DsMC.Branch("MCDstarCt",&MCDstarCt, "MCDstarCt/D");
	
	//----------------------------------------
	//For D* Contamination
	if (debug)cout << "debug D* MC Variables 3--------------------"<< endl;
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
	const int nbinsPt_left = 9;
	double binsPt_left[nbinsPt_left+1] = {4., 5., 6., 7., 8., 12., 16., 24., 40., 100.};
	const int nbinsEta_left = 10;
	double binsEta_left[nbinsEta_left+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1};

	TH1 *hRecPt_MC = makeTH1Rec("hRecPt_MC", nbinsPt_left, binsPt_left );
	TH1 *hRecEta_MC = makeTH1Rec("hRecEta_MC", nbinsEta_left, binsEta_left );
	
	TH1 *hRecPt_Rec1 = makeTH1Rec("hRecPt_Rec1", nbinsPt_left, binsPt_left );
	TH1 *hRecEta_Rec1 = makeTH1Rec("hRecEta_Rec1", nbinsEta_left, binsEta_left );

	TH1 *hRecPt_Rec2 = makeTH1Rec("hRecPt_Rec2", nbinsPt_left, binsPt_left );
	TH1 *hRecEta_Rec2 = makeTH1Rec("hRecEta_Rec2", nbinsEta_left, binsEta_left );

	//----------------------------------------
	// D* REC x GEN
	TH2* hRECxGENDstarMass = makeTH2GenRec( "hRECxGENDstarMass", "Mass[GeV]", 100, 1.94 , 2.1 );
	TH2* hRECxGENDstarPt = makeTH2GenRec( "hRECxGENDstarPt", "pT[GeV]", 100, -1, 40. );
	TH2* hRECxGENDstarEta = makeTH2GenRec( "hRECxGENDstarEta", "#eta", 100, -1, 3. );
	TH2* hRECxGENDstarTime = makeTH2GenRec( "hRECxGENDstarTime", "t[s]", 100, 0, 1*pow(10,-12) );
	TH2* hRECxGENDstarDisXY = makeTH2GenRec( "hRECxGENDstarDisXY", "L[cm]", 100, 0., 0.4 );

	TH2* hRECxGENDstarMass2 = makeTH2GenRec( "hRECxGENDstarMass2", "Mass[GeV]", 100, 1.94 , 2.1 );
	TH2* hRECxGENDstarPt2 = makeTH2GenRec( "hRECxGENDstarPt2", "pT[GeV]", 100, -1, 40. );
	TH2* hRECxGENDstarEta2 = makeTH2GenRec( "hRECxGENDstarEta2", "#eta", 100, -1, 3. );
	TH2* hRECxGENDstarTime2 = makeTH2GenRec( "hRECxGENDstarTime2", "t[s]", 100, 0, 1*pow(10,-12) );
	TH2* hRECxGENDstarDisXY2 = makeTH2GenRec( "hRECxGENDstarDisXY2", "L[cm]", 100, 0., 0.4 );

	//----------------------------------------
	// MC Matching - D* 
	double var_MCD0lifetimeMatching = 0.;
	t_DsMatching.Branch("MCD0lifetimeMatching",&var_MCD0lifetimeMatching, "var_MCD0lifetimeMatching/D");
	double var_deltaRDsMatching = 0.;
	t_DsMatching.Branch("deltaRDsMatching",&var_deltaRDsMatching, "var_deltaRDsMatching/D");

	TH1 *hDstarDeltaR = makeTH1("hDstarDeltaR", 100, 0, 0.06, "#Delta R of the D* ; #Delta R ; Events" );

	TH1 *hDstarDeltaR2 = makeTH1("hDstarDeltaR2", 100, 0, 0.7, "#Delta R of the D* ; #Delta R ; Events" );

	//----------------------------------------
	//For D0 MC
	//-----------------------------------------
	TString MCD0String [] = {"MCpromptD0eta", "MCpromptD0phi", "MCpromptD0pt", "MCpromptD0energy", "MCpromptD0p", "MCpromptD0et", "MCpromptD0rapidity", "MCpromptD0mass", "MCpromptD0_DispAngle", "MCpromptD0lifetime", "MCpromptD0dispXY", "MCpromptD0_Keta", "MCpromptD0_Kpt", "MCpromptD0_Pieta", "MCpromptD0_Pipt" };

	vector< vector<double>* > MCD0Vec;
	vector<double>* MCpromptD0eta =0; vector<double>* MCpromptD0phi =0; vector<double>* MCpromptD0pt =0; vector<double>* MCpromptD0energy =0; vector<double>* MCpromptD0p =0; vector<double>* MCpromptD0et =0; vector<double>* MCpromptD0rapidity =0; vector<double>* MCpromptD0mass =0; vector<double>* MCpromptD0_DispAngle =0; vector<double>* MCpromptD0lifetime =0; vector<double>* MCpromptD0dispXY =0; vector<double>* MCpromptD0_Keta =0; vector<double>* MCpromptD0_Kpt =0; vector<double>* MCpromptD0_Pieta =0; vector<double>* MCpromptD0_Pipt =0; 

	MCD0Vec.push_back(MCpromptD0eta); MCD0Vec.push_back(MCpromptD0phi); MCD0Vec.push_back(MCpromptD0pt); MCD0Vec.push_back(MCpromptD0energy); MCD0Vec.push_back(MCpromptD0p); MCD0Vec.push_back(MCpromptD0et); MCD0Vec.push_back(MCpromptD0rapidity); MCD0Vec.push_back(MCpromptD0mass); MCD0Vec.push_back(MCpromptD0_DispAngle); MCD0Vec.push_back(MCpromptD0lifetime); MCD0Vec.push_back(MCpromptD0dispXY); MCD0Vec.push_back(MCpromptD0_Keta); MCD0Vec.push_back(MCpromptD0_Kpt); MCD0Vec.push_back(MCpromptD0_Pieta); MCD0Vec.push_back(MCpromptD0_Pipt); 

	// Input Branches
	vector< TBranch* > MCD0Branch; MCD0Branch.clear();
	TBranch *b_MCpromptD0eta; TBranch *b_MCpromptD0phi; TBranch *b_MCpromptD0pt; TBranch *b_MCpromptD0energy; TBranch *b_MCpromptD0p; TBranch *b_MCpromptD0et; TBranch *b_MCpromptD0rapidity; TBranch *b_MCpromptD0mass; TBranch *b_MCpromptD0_DispAngle; TBranch *b_MCpromptD0lifetime; TBranch *b_MCpromptD0dispXY; TBranch *b_MCpromptD0_Keta; TBranch *b_MCpromptD0_Kpt; TBranch *b_MCpromptD0_Pieta; TBranch *b_MCpromptD0_Pipt; 

	MCD0Branch.push_back(b_MCpromptD0eta); MCD0Branch.push_back(b_MCpromptD0phi); MCD0Branch.push_back(b_MCpromptD0pt); MCD0Branch.push_back(b_MCpromptD0energy); MCD0Branch.push_back(b_MCpromptD0p); MCD0Branch.push_back(b_MCpromptD0et); MCD0Branch.push_back(b_MCpromptD0rapidity); MCD0Branch.push_back(b_MCpromptD0mass); MCD0Branch.push_back(b_MCpromptD0_DispAngle); MCD0Branch.push_back(b_MCpromptD0lifetime); MCD0Branch.push_back(b_MCpromptD0dispXY); MCD0Branch.push_back(b_MCpromptD0_Keta); MCD0Branch.push_back(b_MCpromptD0_Kpt); MCD0Branch.push_back(b_MCpromptD0_Pieta); MCD0Branch.push_back(b_MCpromptD0_Pipt); 

	//-----------------------------------------
	//ADDRESSING
	for( unsigned int k = 0; k < MCD0Branch.size(); k++)
	{ t1->SetBranchAddress( MCD0String[k], &MCD0Vec[k], &MCD0Branch[k]); }
	//-----------------------------------------
	//Variables 
	double var_MCpromptD0eta =0.; double var_MCpromptD0phi =0.; double var_MCpromptD0pt =0.; double var_MCpromptD0energy =0.; double var_MCpromptD0p =0.; double var_MCpromptD0et =0.; double var_MCpromptD0rapidity =0.; double var_MCpromptD0mass =0.; double var_MCpromptD0_DispAngle =0.; double var_MCpromptD0lifetime =0.; double var_MCpromptD0dispXY =0.; double var_MCpromptD0_Keta =0.; double var_MCpromptD0_Kpt =0.; double var_MCpromptD0_Pieta =0.; double var_MCpromptD0_Pipt =0.;  
	double MCD0Var[] = {var_MCpromptD0eta, var_MCpromptD0phi, var_MCpromptD0pt, var_MCpromptD0energy, var_MCpromptD0p, var_MCpromptD0et, var_MCpromptD0rapidity, var_MCpromptD0mass, var_MCpromptD0_DispAngle, var_MCpromptD0lifetime, var_MCpromptD0dispXY, var_MCpromptD0_Keta, var_MCpromptD0_Kpt, var_MCpromptD0_Pieta, var_MCpromptD0_Pipt };

	sizeArray = sizeof(MCD0Var)/sizeof(MCD0Var[0]); //Number of elements
	vector<double> MCD0VarVec; MCD0VarVec.clear();
	MCD0VarVec.insert(MCD0VarVec.begin(), MCD0Var, MCD0Var+sizeArray );
	//-----------------------------------------'
	// OUTPUT BRANCHES
	for( unsigned int k = 0; k < MCD0VarVec.size(); k++)
	{ t_D0MC.Branch( MCD0String[k], &MCD0VarVec[k], "var_"+MCD0String[k]+"/D"); }
	t_D0MC.Branch("PUWeight",&PUWeight, "PUWeight/D");

	double MCD0Ct =0.;
	t_D0MC.Branch("MCD0Ct",&MCD0Ct, "MCD0Ct/D");

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
	TH1 *hRecD0Pt_MC = makeTH1Rec("hRecD0Pt_MC", nbinsPt_left, binsPt_left );
	TH1 *hRecD0Eta_MC = makeTH1Rec("hRecD0Eta_MC", nbinsEta_left, binsEta_left );
	TH1 *hRecD0Pt_Rec1 = makeTH1Rec("hRecD0Pt_Rec1", nbinsPt_left, binsPt_left );
	TH1 *hRecD0Eta_Rec1 = makeTH1Rec("hRecD0Eta_Rec1", nbinsEta_left, binsEta_left );

	//----------------------------------------
	// D0 REC x GEN
	TH2* hRECxGEND0Mass = makeTH2GenRec( "hRECxGEND0Mass", "Mass[GeV]", 100, 1.76, 1.96 );
	TH2* hRECxGEND0Pt = makeTH2GenRec( "hRECxGEND0Pt", "pT[GeV]", 100, -1, 30. );
	TH2* hRECxGEND0Eta = makeTH2GenRec( "hRECxGEND0Eta", "#eta", 100, -1, 3. );
	TH2* hRECxGEND0Time = makeTH2GenRec( "hRECxGEND0Time", "t[s]", 100, 0, 3*pow(10,-12) );
	TH2* hRECxGEND0DisXY = makeTH2GenRec( "hRECxGEND0DisXY", "L[cm]", 100, 0., 0.4 );

	//----------------------------------------
	// MC Matching - D0 
	double var_MCpromptD0lifetimeMatching = 0.;
	t_D0Matching.Branch("MCpromptD0lifetimeMatching",&var_MCpromptD0lifetimeMatching, "var_MCpromptD0lifetimeMatching/D");
	double var_deltaRD0Matching = 0.;
	t_D0Matching.Branch("deltaRD0Matching",&var_deltaRD0Matching, "var_deltaRD0Matching/D");

	TH1 *hD0DeltaR = makeTH1("hD0DeltaR", 100, 0, 0.06, "#Delta R of the D0 ; #Delta R ; Events" );

#endif

	//----------------------------------------
	// LOOP TREE ENTRIES FOR FILE f1
	for (Long64_t jentry=0; jentry < nentries; jentry++) 
	{
		if (debug)cout << "-------------------- File loop "<< endl;
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
		if (debug)cout << "-------------------- D* GetEntry "<< endl;
		for( unsigned int k = 0; k < DstarVarVec.size(); k++)
		{ DstarBranch[k]->GetEntry(ientry); }
		//----------------------------------------
		//For D* wrong combination
		if (debug)cout << "-------------------- D* Wrong GetEntry "<< endl;
		for( unsigned int k = 0; k < DstarVarVecWrong.size(); k++)
		{ DstarBranchWrong[k]->GetEntry(ientry); }
		//----------------------------------------
		//For D0
		if (debug)cout << "-------------------- D0 GetEntry "<< endl;
		for( unsigned int k = 0; k < D0VarVec.size(); k++)
		{ D0Branch[k]->GetEntry(ientry); }
		//----------------------------------------
		//For D0 Wrong Combination
		if (debug)cout << "-------------------- D0 Wrong GetEntry "<< endl;
		for( unsigned int k = 0; k < D0VarVecWrong.size(); k++)
		{ D0BranchWrong[k]->GetEntry(ientry); }

#if GENvaraibles == 0	
#elif GENvaraibles == 1
		b_PUWeight->GetEntry(ientry);
		//cout << "PUWeight: " << PUWeight << endl;
		//----------------------------------------
		//For D* MC
		if (debug)cout << "-------------------- D* MC GetEntry "<< endl;
		for( unsigned int k = 0; k < MCDstarVarVec.size(); k++)
		{ MCDstarBranch[k]->GetEntry(ientry); }
		b_FlagDstarfromB->GetEntry(ientry);
		//----------------------------------------
		//For D0 MC
		if (debug)cout << "-------------------- D0 MC GetEntry "<< endl;
		for( unsigned int k = 0; k < MCD0VarVec.size(); k++)
		{ MCD0Branch[k]->GetEntry(ientry); }
		b_FlagD0fromB->GetEntry(ientry);
#endif
		
		//----------------------------------------
		// SELECTION
		//----------------------------------------
		//For D* and D0(from D*)
		if (debug)cout << "-------------------- debug For D* and D0(from D*) "<< endl;
		if ( (DstarVec[0]->size() > 0) && FlagDs){
		for(unsigned int i=0; i < DstarVec[0]->size(); i++)
		{  
			//[0]"Dsmass", [1]"Dseta", [2]"Dsphi", [3]"Dspt", [4]"DSDeltaR", [5]"D0mass", [6]"D0eta", [7]"D0phi", [8]"D0pt", [9]"D0fromDSdXY", [10]"D0fromDSsXY", [11]"D0fromDSs3D", [12]"D0fromDSd3D", [13]"Anglephi", [14]"D0_VtxProb", [15]"Dslifetime"
			//PossibleDs++;					
			if( DstarVec[14]->at(i) < 0.01) continue; //DsD0ProbVtx++;
			if( DstarVec[13]->at(i) < 0.99 ) continue; //DsD0Angle++;	
			if( fabs(DstarVec[5]->at(i) - 1.86484) > 0.1 ) continue; //D0fromDsMinusPDG++;
			if( DstarVec[8]->at(i) < 3. ) continue; //DsD0pt++;
			if( ( DstarVec[0]->at(i) - DstarVec[5]->at(i) ) > 0.16) continue; //DsD0deltaM++;

			//Evolution Significance Cut
			double binDstarMass[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
			sizeArray = sizeof(binDstarMass)/sizeof(binDstarMass[0]); // Number of elements
			for( int k = 0; k < sizeArray; k++)
			{ if( DstarVec[11]->at(i) > binDstarMass[k])
			  	{ hDstarMass[k]->Fill(DstarVec[0]->at(i));
				  hDstarLifeTime[k]->Fill(DstarVec[15]->at(i));	}	}


	
			if( DstarVec[11]->at(i) < SigCutD0fromDs) continue; //DsD0Sig++;

			//Evolution pT Cut
			double binDstarMassPt[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
			sizeArray = sizeof(binDstarMassPt)/sizeof(binDstarMassPt[0]); // Number of elements
			for( int k = 0; k < sizeArray; k++)
			{ if( DstarVec[3]->at(i) > binDstarMassPt[k]) { hDstarMassPt[k]->Fill(DstarVec[0]->at(i));}	}

			//----------------------------------------------------
			// Fill Branches
			for( unsigned int k = 0; k < DstarVarVec.size(); k++)
			{	DstarVarVec[k] = DstarVec[k]->at(i); }

			DstarCt = DstarVec[5]->at(i) * DstarVec[9]->at(i) / DstarVec[8]->at(i);
			
			//----------------------------------------------------
			// Fill Histograms
#if GENvaraibles == 0
			for( unsigned int k = 0; k < DstarVarVec.size(); k++)
			{	HistogramsDstar[k]->Fill( DstarVec[k]->at(i) ); }	
#elif GENvaraibles == 1
			for( unsigned int k = 0; k < DstarVarVec.size(); k++)
			{	HistogramsDstar[k]->Fill( DstarVec[k]->at(i), PUWeight ); }
#endif			

		//cout << "Dsmass: "<<  DstarVec[0]->at(i) << " # D0mass: "<<  DstarVec[5]->at(i) << " # Prob: "<<  DstarVec[14]->at(i) << " # AnglePhi: "<<  DstarVec[13]->at(i) << endl;	
			t_analysis.Fill();

			if (debug)cout << "debug D* 12 --------------------"<< endl;
	  	}}

		//----------------------------------------
		//For D* Wrong combination
		if (debug)cout << "-------------------- debug For D* Wrong "<< endl;
		if ( (DstarVecWrong[0]->size() > 0) && FlagDs){
		for(unsigned int i=0; i < DstarVecWrong[0]->size(); i++)
		{  		
			//[0]"DsmassWrong", [1]"DsetaWrong", [2]"DsphiWrong", [3]"DsptWrong", [4]"D0massWrong", [5]"D0etaWrong", [6]"D0phiWrong", [7]"D0ptWrong", [8]"D0fromDSdXYWrong", [9]"D0fromDSsXYWrong", [10]"AnglephiWrong", [11]"D0_VtxProbWrong", [12]"DslifetimeWrong"
			//PossibleDs++;

			if( DstarVecWrong[11]->at(i) < 0.01) continue; //DsD0ProbVtx++;		
			if( DstarVecWrong[10]->at(i) < 0.99 ) continue; //DsD0Angle++;	
			if( fabs(DstarVecWrong[4]->at(i) - 1.86484) > 0.1 ) continue; //D0fromDsMinusPDG++;
			if( DstarVecWrong[3]->at(i) < 3. ) continue; //DsD0pt++;
			if( (DstarVecWrong[0]->at(i) - DstarVecWrong[4]->at(i)) > 0.16) continue; //DsD0deltaM++;*/
			if( DstarVecWrong[9]->at(i) < SigCutD0fromDs) continue; //DsD0Sig++;

			//----------------------------------------------------
			// Fill Branches
			for( unsigned int k = 0; k < DstarVecWrong.size(); k++)
			{	DstarVarVecWrong[k] = DstarVecWrong[k]->at(i); }

			DstarCtWrong = DstarVecWrong[4]->at(i) * DstarVecWrong[8]->at(i) / DstarVecWrong[7]->at(i);

			//----------------------------------------------------
			// Fill Histograms
#if GENvaraibles == 0
			for( unsigned int k = 0; k < DstarVecWrong.size(); k++)
			{	HistogramsDstarWrong[k]->Fill( DstarVecWrong[k]->at(i) ); }	
#elif GENvaraibles == 1
			for( unsigned int k = 0; k < DstarVecWrong.size(); k++)
			{	HistogramsDstarWrong[k]->Fill( DstarVecWrong[k]->at(i), PUWeight ); }
#endif

			t_DstarWrongCombination.Fill();
			if (debug)cout << "-------------------- End For D* Wrong "<< endl;
	  	}}
		//----------------------------------------
		//For D0 
		//-----------------------------------------
		if (debug)cout << "-------------------- debug For D0 "<< endl;
		if ( (D0Vec[0]->size() > 0) && FlagD0){
		for(unsigned int i=0; i < D0Vec[0]->size(); i++)
		{	
			//[0]"D0Kpimass", [1]"D0Kpipt", [2]"D0Kpieta", [3]"D0Kpiphi", [4]"D0lifetime", [5]"D0Kpi_VtxProb", [6]"D0KpiDispAngle", [7]"D0KpisXY", [8]"D0KpidXY", [9]"D0Kpis3D", [10]"D0Kpid3D"			
			//PossibleD0++;
			if( D0Vec[5]->at(i) < 0.01) continue; //D0ProbVtx++;	
			if( D0Vec[6]->at(i) < 0.99 ) continue; //D0Angle++;
			if( abs(D0Vec[0]->at(i)-1.86484) > 0.1) continue;	//CountD0minusPDG++;

			//Evolution pT Cut
			if( D0Vec[9]->at(i) < SigCutD0)
			{	double binD0Pt[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
				sizeArray = sizeof(binD0Pt)/sizeof(binD0Pt[0]); // Number of elements
				for( int k = 0; k < sizeArray; k++)
				{ if( D0Vec[1]->at(i)> binD0Pt[k]) { hD0MassPt[k]->Fill(D0Vec[0]->at(i));}	}	
			}

			if( D0Vec[1]->at(i) < 2. ) continue; //CountD0pt++;

			//Evolution Significance Cut
			double binD0Sig[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
			sizeArray = sizeof(binD0Sig)/sizeof(binD0Sig[0]); // Number of elements
			for( int k = 0; k < sizeArray; k++)
			{ if( D0Vec[9]->at(i)>binD0Sig[k]) { hD0Mass[k]->Fill(D0Vec[0]->at(i));}	}	
			
			if( D0Vec[9]->at(i) < SigCutD0) continue; //D0Sig++;

			//----------------------------------------------------
			// Fill Branches
			for( unsigned int k = 0; k < D0VarVec.size(); k++)
			{ D0VarVec[k] = D0Vec[k]->at(i); }

			D0Ct = D0Vec[0]->at(i)*D0Vec[8]->at(i)/D0Vec[1]->at(i);
			
			//----------------------------------------------------
			// Fill Histograms
#if GENvaraibles == 0
			for( unsigned int k = 0; k < D0VarVec.size(); k++)
			{	HistogramsD0[k]->Fill( D0Vec[k]->at(i) ); }	
#elif GENvaraibles == 1
			for( unsigned int k = 0; k < D0VarVec.size(); k++)
			{	HistogramsD0[k]->Fill( D0Vec[k]->at(i), PUWeight ); }
#endif	
			t_D0analysis.Fill();
			if (debug)cout << "-------------------- End For D0 "<< endl;
		}}

		//----------------------------------------
		//For D0 Wrong Combination
		//-----------------------------------------
		if (debug)cout << "-------------------- debug For D0 "<< endl;
		if ( (D0VecWrong[0]->size() > 0) && FlagD0){
		for(unsigned int i=0; i < D0VecWrong[0]->size(); i++)
		{	
			//[0]"D0KpimassWrong", [1]"D0KpiptWrong", [2]"D0KpietaWrong", [3]"D0KpiphiWrong", [4]"D0lifetimeWrong", [5]"D0Kpi_VtxProbWrong", [6]"D0KpiDispAngleWrong", [7]"D0KpisXYWrong", [8]"D0KpidXYWrong", [9]"D0Kpis3DWrong", [10]"D0Kpid3DWrong", [11]"D0KpikTWrong",		
			//PossibleD0++;
			if( D0VecWrong[5]->at(i) < 0.01) continue; //D0ProbVtx++;	
			if( D0VecWrong[6]->at(i) < 0.99 ) continue; //D0Angle++;
			if( abs(D0VecWrong[0]->at(i)-1.86484) > 0.1) continue;	//CountD0minusPDG++;
			if( D0VecWrong[1]->at(i) < 2. ) continue; //CountD0pt++;		
			if( D0VecWrong[9]->at(i) < SigCutD0) continue; //D0Sig++;

			//----------------------------------------------------
			// Fill Branches
			for( unsigned int k = 0; k < D0VarVecWrong.size(); k++)
			{ D0VarVecWrong[k] = D0VecWrong[k]->at(i); }

			D0CtWrong = D0VecWrong[0]->at(i)*D0VecWrong[8]->at(i)/D0VecWrong[1]->at(i);
			
			//----------------------------------------------------
			// Fill Histograms
#if GENvaraibles == 0
			for( unsigned int k = 0; k < D0VarVecWrong.size(); k++)
			{	HistogramsD0Wrong[k]->Fill( D0VecWrong[k]->at(i) ); }	
#elif GENvaraibles == 1
			for( unsigned int k = 0; k < D0VarVecWrong.size(); k++)
			{	HistogramsD0Wrong[k]->Fill( D0VecWrong[k]->at(i), PUWeight ); }
#endif	
			t_D0analysis.Fill();
			if (debug)cout << "-------------------- End For D0 "<< endl;
		}}

#if GENvaraibles == 0	
#elif GENvaraibles == 1
		//----------------------------------------
		//For D* MC
		if (debug)cout << "-------------------- debug For D* MC "<< endl;
		if ( (MCDstarVec[0]->size() > 0) && FlagDsMC){ //D* MC protection
		for(unsigned int i=0; i < MCDstarVec[0]->size(); i++) // loop D* MC
		{
			//[0]"MCDseta", [1]"MCDsphi", [2]"MCDspt", [3]"MCDsenergy", [4]"MCDsp", [5]"MCDset", [6]"MCDsmass", [7]"MCD0eta", [8]"MCD0phi", [9]"MCD0pt", [10]"MCD0energy", [11]"MCD0p", [12]"MCD0et", [13]"MCD0rapidity", [14]"MCD0mass", [15]"MCD0lifetime", [16]"MCD0dispXY", [17]"MCDsSeta", [18]"MCDsSpt", [19]"MCDsKeta", [20]"MCDsKpt", [21]"MCDsPieta", [22]"MCDsPipt", 

			// Same cinematic region of Reconstructed tracks 
			if ( fabs( MCDstarVec[17]->at(i) ) > 2.4 ) continue; //SlowPion
			if ( fabs( MCDstarVec[19]->at(i) ) > 2.4 ) continue; //Kaon	
			if ( fabs( MCDstarVec[21]->at(i) ) > 2.4 ) continue; //Pion
			if ( MCDstarVec[18]->at(i) < 0.5 ) continue; //SlowPion
			if ( MCDstarVec[20]->at(i) < 0.5 ) continue; //Kaon
			if ( MCDstarVec[22]->at(i) < 0.5 ) continue; //Pion
	
			// Cinematic Region For DO
			//if ( MCDstarVec[9]->at(i) < 3. ) continue;
			// Cinematic Region For Eta and pT of D* and DO
			if ( fabs(MCDstarVec[0]->at(i)) > 2.1 ) continue;
			if (  (MCDstarVec[2]->at(i) < 4. ) or (MCDstarVec[2]->at(i) > 100.) ) continue;

			// Fill Branches
			for( unsigned int k = 0; k < MCDstarVarVec.size(); k++)
			{ MCDstarVarVec[k] = MCDstarVec[k]->at(i); }

			//Proper Decay Time
			MCDstarCt = MCDstarVec[6]->at(i) * MCDstarVec[16]->at(i) / MCDstarVec[9]->at(i);

			// D* recosntruction Efficiency pT and Eta
			hRecEta_MC->Fill( abs(MCDstarVec[0]->at(i)) );
			hRecPt_MC->Fill( MCDstarVec[2]->at(i) );

			t_DsMC.Fill();
		
			//----------------------------------------
			//For D* Matching
			//----------------------------------------
			if (debug)cout << "-------------------- debug For D* Matching "<< endl;
			if ( DstarVec[0]->size() > 0 ){
			double minDeltaR = 1000.;
			int id = -1; // Protection
			double deltaR = 0.;
			double deltaEta= 0.;
			double deltaPhi= 0.;
			for(unsigned int j=0; j < DstarVec[0]->size(); j++)
			{	
			//[0]"Dsmass", [1]"Dseta", [2]"Dsphi", [3]"Dspt", [4]"DSDeltaR", [5]"D0mass", [6]"D0eta", [7]"D0phi", [8]"D0pt", [9]"D0fromDSdXY", [10]"D0fromDSsXY", [11]"D0fromDSs3D", [12]"D0fromDSd3D", [13]"Anglephi", [14]"D0_VtxProb", [15]"Dslifetime"
				if( DstarVec[14]->at(j) < 0.01) continue;	
				if( DstarVec[13]->at(j) < 0.99 ) continue; 
				if( fabs(DstarVec[5]->at(j) - 1.86484) > 0.1 ) continue;
				if( DstarVec[8]->at(j) < 3. ) continue;
				if( (DstarVec[0]->at(j) - DstarVec[5]->at(j)) > 0.16) continue;
				if( DstarVec[11]->at(j) < SigCutD0fromDs) continue; //DsD0Sig++;
						
				//----------------------------------------
				//For D* Contamination
				if (debug)cout << "debug For D* Contamination --------------------"<< endl;
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

				var_MCD0lifetimeMatching = MCDstarVec[15]->at(i);
				var_deltaRDsMatching = deltaR;
				hDstarDeltaR->Fill( deltaR );

				hDstarDeltaR2->Fill( deltaR );
				
				double recMass = -1; double GenMass = -1; // Protection
				double recPt = -1; double GenPt = -1;
				double recEta = -1; double GenEta = -1;
				double recTime = -1; double GenTime = -1;
				double recDisXY = -1; double GenDisXY = -1;
				//Histograms for efficiency
				if (debug)cout << "debug For D* efficiency --------------------"<< endl;
				if( deltaR < 0.03){

					// D* recosntruction Efficiency pT and Eta
					hRecEta_Rec1->Fill( abs(DstarVec[1]->at(id)) );
					hRecPt_Rec1->Fill( DstarVec[3]->at(id));
		
					// D* REC x GEN
					recMass = DstarVec[0]->at(id); GenMass = MCDstarVec[6]->at(i);
					recPt = DstarVec[3]->at(id); GenPt = MCDstarVec[2]->at(i);
					recEta = DstarVec[1]->at(id); GenEta = MCDstarVec[0]->at(i);
					recTime = DstarVec[15]->at(id); GenTime = MCDstarVec[15]->at(i);
					recDisXY = DstarVec[9]->at(id); GenDisXY = MCDstarVec[16]->at(i);	
					hRECxGENDstarMass->Fill( recMass, GenMass );
					hRECxGENDstarPt->Fill( recPt, GenPt );
					hRECxGENDstarEta->Fill( abs(recEta), abs(GenEta) );
					hRECxGENDstarTime->Fill( recTime, GenTime );
					hRECxGENDstarDisXY->Fill( recDisXY, GenDisXY );					
				}
				if( deltaR < 0.7){

					// D* recosntruction Efficiency pT and Eta
					hRecEta_Rec2->Fill( abs(DstarVec[1]->at(id)) );
					hRecPt_Rec2->Fill( DstarVec[3]->at(id));
			
					// D* REC x GEN
					recMass = DstarVec[0]->at(id); GenMass = MCDstarVec[6]->at(i);
					recPt = DstarVec[3]->at(id); GenPt = MCDstarVec[2]->at(i);
					recEta = DstarVec[1]->at(id); GenEta = MCDstarVec[0]->at(i);
					recTime = DstarVec[15]->at(id); GenTime = MCDstarVec[15]->at(i);
					recDisXY = DstarVec[9]->at(id); GenDisXY = MCDstarVec[16]->at(i);	
					hRECxGENDstarMass2->Fill( recMass, GenMass );
					hRECxGENDstarPt2->Fill( recPt, GenPt );
					hRECxGENDstarEta2->Fill( abs(recEta), abs(GenEta) );
					hRECxGENDstarTime2->Fill( recTime, GenTime );
					hRECxGENDstarDisXY2->Fill( recDisXY, GenDisXY );					
				}
				t_DsMatching.Fill();	
			}
			} // End D* Matching

		} // End loop D* MC
		} // End D* MC protection
		//----------------------------------------
		//For D0 MC
		if (debug)cout << "debug For D0 MC --------------------"<< endl;
		if ( (MCD0Vec[0]->size() > 0) && FlagD0MC){ //D0 MC protection
		for(unsigned int i=0; i < MCD0Vec[0]->size(); i++) // loop D0 MC
		{
			//[0]"MCpromptD0eta", [1]"MCpromptD0phi", [2]"MCpromptD0pt", [3]"MCpromptD0energy", [4]"MCpromptD0p", [5]"MCpromptD0et", [6]"MCpromptD0rapidity", [7]"MCpromptD0mass", [8]"MCpromptD0_DispAngle", [9]"MCpromptD0lifetime", [10]"MCpromptD0dispXY", [11]"MCpromptD0_Keta", [12]"MCpromptD0_Kpt", [13]"MCpromptD0_Pieta", [14]"MCpromptD0_Pipt"

			// Same cinematic region of tracks ( Reconstruction Cuts )
			if ( fabs( MCD0Vec[11]->at(i)) > 2.4 ) continue;
			if ( fabs( MCD0Vec[12]->at(i)) > 2.4 ) continue;
			if ( MCD0Vec[13]->at(i) < 0.8 ) continue;
			if ( MCD0Vec[14]->at(i) < 0.8 ) continue;
			
			// Cinematic Region For Eta and pT of D0
			if ( fabs(MCD0Vec[0]->at(i)) > 2.1 ) continue;
			if (  (MCD0Vec[2]->at(i) < 4 ) or (MCD0Vec[2]->at(i) > 100) ) continue;		

			// Fill Branches
			for( unsigned int k = 0; k < MCD0VarVec.size(); k++)
			{ MCD0VarVec[k] = MCD0Vec[k]->at(i); }

			MCD0Ct = MCD0Vec[7]->at(i) * MCD0Vec[10]->at(i) / MCD0Vec[2]->at(i);

			// D0 recosntruction Efficiency pT and Eta
			hRecD0Eta_MC->Fill(abs(MCD0Vec[0]->at(i)));
			hRecD0Pt_MC->Fill(MCD0Vec[2]->at(i));

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
				//[0]"D0Kpimass", [1]"D0Kpipt", [2]"D0Kpieta", [3]"D0Kpiphi", [4]"D0lifetime", [5]"D0Kpi_VtxProb", [6]"D0KpiDispAngle", [7]"D0KpisXY", [8]"D0KpidXY", [9]"D0Kpis3D", [10]"D0Kpid3D"
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
				//{"MCpromptD0eta", "MCpromptD0phi", "MCpromptD0pt", "MCpromptD0energy", "MCpromptD0p", "MCpromptD0et", "MCpromptD0rapidity", "MCpromptD0mass", "MCpromptD0_DispAngle", "MCpromptD0lifetime", "MCpromptD0dispXY" }
				deltaEta = MCD0Vec[0]->at(i) - D0Vec[2]->at(j);
				deltaPhi = MCD0Vec[1]->at(i) - D0Vec[3]->at(j);
				deltaR = sqrt( pow( deltaEta, 2) + pow( deltaPhi ,2));
				if( deltaR < minDeltaR) { id = j; minDeltaR = deltaR;} //Get the lower delta
				if (debug)cout << "debug D* 20 --------------------"<< endl;
			}

			if (id != -1){

				var_MCpromptD0lifetimeMatching = MCD0Vec[9]->at(i);		
				var_deltaRD0Matching = deltaR;
				hD0DeltaR->Fill( deltaR );

				double recMass = -1; double GenMass = -1; // Protection
				double recPt = -1; double GenPt = -1;
				double recEta = -1; double GenEta = -1;
				double recTime = -1; double GenTime = -1;
				double recDisXY = -1; double GenDisXY = -1;
				//Histograms for efficiency
				if( deltaR < 0.03){
							
					// D* recosntruction Efficiency pT and Eta					
					hRecD0Pt_Rec1->Fill( D0Vec[1]->at(id) );
					hRecD0Eta_Rec1->Fill( abs(D0Vec[2]->at(id)) );

					recMass = D0Vec[0]->at(id); GenMass = MCD0Vec[7]->at(i);
					recPt = D0Vec[1]->at(id); GenPt = MCD0Vec[2]->at(i);
					recEta = D0Vec[2]->at(id);	GenEta = MCD0Vec[0]->at(i);
					recTime = D0Vec[4]->at(id);	GenTime = MCD0Vec[9]->at(i);
					recDisXY = D0Vec[8]->at(id); GenDisXY = MCD0Vec[10]->at(i);

					hRECxGEND0Mass->Fill( recMass, GenMass );
					hRECxGEND0Pt->Fill( recPt, GenPt );
					hRECxGEND0Eta->Fill( abs(recEta), abs(GenEta) );
					hRECxGEND0Time->Fill( recTime, GenTime );
					hRECxGEND0DisXY->Fill( recDisXY, GenDisXY );	
				}
				
				t_D0Matching.Fill();	
			} // End ID
			} // End D* Matching


		} // End loop D0 MC
		} // End D0 MC protection
		
#endif

	if (debug)cout << "debug 24 --------------------"<< endl;	
	}//End loop tree entries for file f1


	TLatex* tex13;
	TCanvas* canvas = new TCanvas("PreliminarStudyDs","",900,600);

#if GENvaraibles == 0
#elif GENvaraibles == 1
	TString titleRec;
	// Rec D* Efficiency
	titleRec = "Reconstruction Efficiency D*; pT; Efficiency" ;
	TGraphAsymmErrorsRec( hRecPt_Rec1, hRecPt_MC, "hRecPt_Rec1", titleRec, nbinsPt_left, Dataset, path2 );
	titleRec = "Reconstruction Efficiency D*; #eta; Efficiency" ;
	TGraphAsymmErrorsRec( hRecEta_Rec1, hRecEta_MC, "hRecEta_Rec1", titleRec, nbinsEta_left, Dataset, path2 );

	titleRec = "Reconstruction Efficiency D*; pT; Efficiency" ;
	TGraphAsymmErrorsRec( hRecPt_Rec2, hRecPt_MC, "hRecPt_Rec2", titleRec, nbinsPt_left, Dataset, path2 );
	titleRec = "Reconstruction Efficiency D*; #eta; Efficiency" ;
	TGraphAsymmErrorsRec( hRecEta_Rec2, hRecEta_MC, "hRecEta_Rec2", titleRec, nbinsEta_left, Dataset, path2 );

	// Rec D0 Efficiency
	titleRec = "Reconstruction Efficiency D0; pT; Efficiency" ;
	TGraphAsymmErrorsRec( hRecD0Pt_Rec1, hRecD0Pt_MC, "hRecD0Pt_Rec1", titleRec, nbinsPt_left, Dataset, path2 );	
	titleRec = "Reconstruction Efficiency D0; #eta; Efficiency" ;
	TGraphAsymmErrorsRec( hRecD0Eta_Rec1, hRecD0Eta_MC, "hRecD0Eta_Rec1", titleRec, nbinsEta_left, Dataset, path2 );

	//**********************************************************
	canvas = new TCanvas("StudyRecGenDstar","",900,600); canvas->Divide(3,2);
	TH1* nRecGenDstar[] = {hRECxGENDstarMass, hRECxGENDstarPt, hRECxGENDstarEta, hRECxGENDstarTime, hRECxGENDstarDisXY };
	sizeArray = sizeof(nRecGenDstar) / sizeof(nRecGenDstar[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ ){	
	canvas->cd(k+1);	nRecGenDstar[k]->Draw(); tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.03");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw(); }
	canvas->SaveAs(Form("%s%s/StudyRecGenDstar_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/StudyRecGenDstar_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 37 --------------------"<< endl;
	//**********************************************************
	canvas = new TCanvas("StudyRecGenDstar","",900,600); canvas->Divide(3,2);
	TH1* nRecGenDstar2[] = {hRECxGENDstarMass2, hRECxGENDstarPt2, hRECxGENDstarEta2, hRECxGENDstarTime2, hRECxGENDstarDisXY2 };
	sizeArray = sizeof(nRecGenDstar2) / sizeof(nRecGenDstar[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ ){	
	canvas->cd(k+1);	nRecGenDstar[k]->Draw(); tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.03");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw(); }
	canvas->SaveAs(Form("%s%s/StudyRecGenDstar_%sDelta07.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/StudyRecGenDstar_%sDelta07.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 37 --------------------"<< endl;

	//**********************************************************
	canvas = new TCanvas("StudyRecGenD0","",900,600); canvas->Divide(3,2);
	TH1* nRecGenD0[] = {hRECxGEND0Mass, hRECxGEND0Pt, hRECxGEND0Eta, hRECxGEND0Time, hRECxGEND0DisXY };
	sizeArray = sizeof(nRecGenD0) / sizeof(nRecGenD0[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ ){	
	canvas->cd(k+1);	nRecGenD0[k]->Draw(); tex13 = new TLatex(0.15, 0.855, "#Delta R < 0.03");
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw(); }
	canvas->SaveAs(Form("%s%s/StudyRecGenD0_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/StudyRecGenD0_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 38 --------------------"<< endl;


	/*if (Dataset == "MC_DStarToD0Pi_D0KPi" )
	{
		canvas = new TCanvas("CanvasDstarDeltaR","",900,600);
		hDstarDeltaR->Draw("Histo");
		canvas->SaveAs(Form("%s%s/DeltaRDstar.pdf", path2.c_str(), Dataset.c_str() ) );
		canvas->SaveAs(Form("%s%s/DeltaRDstar.png", path2.c_str(), Dataset.c_str() ) );

		canvas = new TCanvas("CanvasDsDeltaR","",900,600);
		hD0DeltaR->Draw("Histo");
		canvas->SaveAs(Form("%s%s/DeltaRD0.pdf", path2.c_str(), Dataset.c_str() ));
		canvas->SaveAs(Form("%s%s/DeltaRD0.png", path2.c_str(), Dataset.c_str() ));
	}*/
	
#endif

	//-----------------------------------------------------------------------------
	canvas = new TCanvas("PreliminarStudyDsMass","",900,600); canvas->Divide(3,3);	
	TString DsmassArray[] = { "= 0", "> 0.5", "> 1.0", ">1.5", "> 2.0", "> 2.5", "> 3.0"};
	sizeArray = sizeof(DsmassArray) / sizeof(DsmassArray[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ )
	{ 	canvas->cd(k+1);	hDstarMass[k]->Draw("e1p");
		tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig "+DsmassArray[k]+") }");
		if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig "+DsmassArray[k]+") }");}
		if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias pp-13Tev (Sig "+DsmassArray[k]+") }");}
		if(Dataset == "MC_BdToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BdToDStarX_ToD0Pi pp-13Tev (Sig "+DsmassArray[k]+") }");}
		if(Dataset == "MC_BuToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BuToDStarX_ToD0Pi pp-13Tev (Sig "+DsmassArray[k]+") }");}
		tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();	
	}
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMass_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyDsMass_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 32 --------------------"<< endl;
//-----------------------------------------------------------------------------
	canvas = new TCanvas("PreliminarStudyDsLifetime","",900,600); canvas->Divide(3,3);	
	for( int k=0; k < sizeArray; k++ )
	{ 	canvas->cd(k+1);	hDstarLifeTime[k]->Draw("e1p");
		tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (Sig "+DsmassArray[k]+") }");
		if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (Sig "+DsmassArray[k]+") }");}
		if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias_ pp-13Tev (Sig "+DsmassArray[k]+") }");}
		if(Dataset == "MC_BdToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BdToDStarX_ToD0Pi pp-13Tev (Sig "+DsmassArray[k]+") }");}
		if(Dataset == "MC_BuToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BuToDStarX_ToD0Pi pp-13Tev (Sig "+DsmassArray[k]+") }");}
		tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();	 
	}
	canvas->SaveAs(Form("%s%s/StudyDsLifetime_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/StudyDsLifetime_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	if (debug)cout << "debug 32 --------------------"<< endl;
	//-----------------------------------------------------------------------------	
	canvas = new TCanvas("PreliminarStudyDsMassPt","",900,600); canvas->Divide(3,3);	
	TString DsmassPtArray[] = { "= 0", "> 0.5", "> 1.0", ">1.5", "> 2.0", "> 2.5", "> 3.0", "> 3.5", "> 4.0"};
	sizeArray = sizeof(DsmassPtArray) / sizeof(DsmassPtArray[0]); // Number of elements
	for( int k=0; k < sizeArray; k++ )
	{ 	canvas->cd(k+1);	hDstarMassPt[k]->Draw("e1p");
		tex13 = new TLatex(0.15, 0.855, "#bf{DATA-Bparking pp-13Tev (pT "+DsmassPtArray[k]+") }");
		if(Dataset == "MC_DStarToD0Pi_D0KPi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-DStar-pythia8-evtgen pp-13Tev (pT "+DsmassPtArray[k]+") }");}
		if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias pp-13Tev (pT "+DsmassPtArray[k]+") }");}
		if(Dataset == "MC_BdToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BdToDStarX_ToD0Pi pp-13Tev (pT "+DsmassPtArray[k]+") }");}
		if(Dataset == "MC_BuToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BuToDStarX_ToD0Pi pp-13Tev (pT "+DsmassPtArray[k]+") }");}
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
		if(Dataset == "MC_MinBias"){tex13 = new TLatex(0.15, 0.855, "#bf{MC-MinBias pp-13Tev (pT "+D0KpimassArray[k]+") }");}
		if(Dataset == "MC_BdToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BdToDStarX_ToD0Pi pp-13Tev (pT "+D0KpimassArray[k]+") }");}
		if(Dataset == "MC_BuToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BuToDStarX_ToD0Pi pp-13Tev (pT "+D0KpimassArray[k]+") }");}
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
		if(Dataset == "MC_BdToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BdToDStarX_ToD0Pi pp-13Tev (pT "+D0KpimassPtArray[k]+") }");}
		if(Dataset == "MC_BuToDStarX_ToD0Pi"){tex13 = new TLatex(0.15, 0.855, "#bf{MC_BuToDStarX_ToD0Pi pp-13Tev (pT "+D0KpimassPtArray[k]+") }");}
		tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04); tex13->Draw();	 }
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0MassPt_%s.pdf", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));
	canvas->SaveAs(Form("%s%s/PreliminarStudyD0MassPt_%s.png", path2.c_str(), Dataset.c_str(), Dataset.c_str() ));		
	if (debug)cout << "debug 35 --------------------"<< endl;
	//-----------------------------------------------------------------------------
	
	//Write in the root file
	f_analysis.cd();

	t_analysis.Write();
	t_DstarWrongCombination.Write();
	t_D0analysis.Write();
#if GENvaraibles == 0	
#elif GENvaraibles == 1
	hDstarDeltaR->Write( );
	hDstarDeltaR2->Write( );
	hD0DeltaR->Write( );
	t_DsMC.Write(); //Write() if a file is open, this function writes a root objectics on it.
	t_D0MC.Write();
	t_DsMatching.Write();
	t_D0Matching.Write();
	t_DstarContamination.Write();
	t_D0Contamination.Write();
#endif	

	// Write Histograms
	sizeArray = sizeof(HistogramsDstar) / sizeof(HistogramsDstar[0]); // Number of elements
	for( int k = 0; k < sizeArray ; k++)
	{	HistogramsDstar[k]->Write( ); }

	sizeArray = sizeof(HistogramsDstarWrong) / sizeof(HistogramsDstarWrong[0]); // Number of elements
	for( int k = 0; k < sizeArray ; k++)
	{	HistogramsDstarWrong[k]->Write( ); }

	sizeArray = sizeof(HistogramsD0) / sizeof(HistogramsD0[0]); // Number of elements
	for( int k = 0; k < sizeArray ; k++)
	{	HistogramsD0[k]->Write( ); }
	
	sizeArray = sizeof(HistogramsD0Wrong) / sizeof(HistogramsD0Wrong[0]); // Number of elements
	for( int k = 0; k < sizeArray ; k++)
	{	HistogramsD0Wrong[k]->Write( ); }	


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
TH1* makeTH1(const char* name, Double_t BIN, Double_t NMIN, Double_t NMAX, const char* TITLES ) 
{
   //Create a TH1D (one dimension) histogram
   TH1D* hh = new TH1D(name, name, BIN, NMIN, NMAX) ;
	hh->SetTitle(TITLES); hh->SetMarkerStyle(7);
	hh->Sumw2();
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
void TGraphAsymmErrorsRec( TH1* h1, TH1* h2, const char* NAME, const char* hTITLE, int Nbins, string DATASETNAME, string path )
{	//Make a TGraphAsymmErrors using histograms as input
	if( TEfficiency::CheckConsistency(*h1,*h2) )
	{	
		TCanvas* canvas = new TCanvas("canvas","",1200,600);
		TGraphAsymmErrors *gr = new TGraphAsymmErrors(h1,h2);
		gr->SetTitle(hTITLE); 	gr->SetMarkerColor(4);
		gr->SetMarkerStyle(21); gr->Draw("AP");
		canvas->SaveAs(Form("%s%s/TGraphAsymmErrors%s.pdf", path.c_str(), DATASETNAME.c_str(), NAME));
		canvas->SaveAs(Form("%s%s/TGraphAsymmErrors%s.png", path.c_str(), DATASETNAME.c_str(), NAME));
		// Divide
		TH1 *h1clone = (TH1F*) h1->Clone();
		h1clone->Divide(h2);
		// Be careful - The first bin is difeten for each histogram *gr and *h1,*h2
		for( int k=0; k < Nbins; k++ ){	
			cout<< " Bin " << k << ": " << h2->GetBinContent(k+1) << " & " << h1->GetBinContent(k+1) << " & "  << h1clone->GetBinContent(k+1) << " pm "<< gr->GetErrorY(k) << endl;
		}
		delete canvas;
	}
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
