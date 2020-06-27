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
	bool debug = false;
		
	//TChain* t1 = callTchain("/eos/user/r/ragomesd/crab/ParkingBPH2/ParkingBPH_Run2018A_MINIAOD/200105_042547/0000/");
	TChain* t1 = callTchain("./");

	//**********************************************************		
	Long64_t nentries = t1->GetEntries(); //Reading Number of tree entries of the file
	nentries = 5000; //Test
	cout<< "Number of tree entries: "<< nentries << endl;
	Long64_t partpercent = nentries*0.05; //Percent done to show
	//double c = 299792458; //Light speed [m/s]
	
	//-------Reading the root file-------------------------------------	
	//TFile *f1 = new TFile("D0DstarDataBparking_7.root");
	//TTree *t1 = (TTree*)f1->Get("analysis/data");
	
	//-------Create root file-------------------------------------
	//TFile f_analysis("DMesonsHistograms.root","recreate"); //Creates root file
	TFile f_analysis(Form("%s%s.root", path.c_str(), Dataset.c_str() ),"recreate"); //Creates root file
	TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree

	//-----------------------------------------------------
	//V A R I A B L E S
	//-----------------------------------------------------
	//It must be equal to "0", not "0.", even it is a double. Possible segmental error if you compile it with c++.

	//----------------------------------------
	//For D*
	//-----------------------------------------
	vector<TString> StringDstar_vec;
	TString StringArray [] = {"Dsmass"};
	int sizeArray = sizeof(StringArray)/sizeof(StringArray[0]); 
	StringDstar_vec.insert(StringDstar_vec.begin(), StringArray, StringArray+sizeArray );

	vector<double> *Dsmass;
	//vector<vector<double>*> *Dstar_vec;
	//Dstar_vec.push_back(Dsmass);

	vector<double> var_Dstar_vec;
	double var_Dsmass;
	var_Dstar_vec.push_back(var_Dsmass);
	for( unsigned int k = 0; k < StringDstar_vec.size(); k++)
	{ t_analysis.Branch( StringDstar_vec[k], &var_Dstar_vec, "var_"+StringDstar_vec[k]+"/D"); } 
	
	TBranch *b_Dsmass; 
	t1->SetBranchAddress( StringDstar_vec[0], &Dsmass, &b_Dsmass);

	
	//TBranch *b_Dsmass; t1->SetBranchAddress("Dsmass",&Dsmass,&b_Dsmass);

	for (Long64_t jentry=0; jentry < nentries; jentry++) //loop tree entries for file f1
	{
		Long64_t ientry = t1->LoadTree(jentry);
      //cout << "nentries " << nentries << " ientry " << ientry << " jentry " << jentry <<endl;
      if (ientry < 0) break;

		//Output about percent program executed
		double percent = (jentry*100)/nentries;		
		if ( jentry % partpercent == 0) cout <<"*****"<< percent << "per cent done***" << endl;
		//cout <<"*****entrada"<< jentry << "***" << endl;
		
		b_Dsmass->GetEntry(ientry);
		
		//----------------------------------------
		//For D* and D0(from D*)
		//-----------------------------------------
		if ( (Dsmass->size() > 0) && FlagDs){
		for(unsigned int i=0; i < Dsmass->size(); i++)
		{  
			var_Dsmass = Dsmass->at(i);
			cout << "Dsmass: "<<  Dsmass->at(i) << endl;
			t_analysis.Fill();

			if (debug)cout << "debug D* 12 --------------------"<< endl;
  		}
		}
		
		
		if (debug)cout << "debug 24 --------------------"<< endl;	
	}//End loop tree entries for file f1

	//t_analysis.Branch("D0mass",&D0mass);
	//TBranch *D0mass_branch; D0mass_branch = 
	f_analysis.cd();

	t_analysis.Write();  //Write in the root file


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
