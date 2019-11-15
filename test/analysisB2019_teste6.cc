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


//gSystem->AddIncludePath("/home/raphael/Desktop/analysis_2018");

using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void analysisB2019_teste6()
{
	//call a file for a histogram style (Optional)
	//gROOT->LoadMacro("styleTDR.C"); 
	//setTDRStyle();

	bool debug = true;
	if (debug)cout << "Teste 1 --------------------"<< endl;
	std::vector<double>* D0mass;
	
	//std::vector<double>* Trkpipt = 0.;
	
	if (debug)cout << "Teste 2 --------------------"<< endl;
	//std::vector<double> vectorInvariantMass_Dstar = 0.;
	//std::vector<double> vectorInvariantMass_D0 = 0.;
	
	//*****Creating Histgrams******************************************************	

	//Histogramas cinematics quantities of the D0
	TH1D *D0mass_Histo = new TH1D("D0mass_Histo","D0mass_Histo",100,1.7,2.0);
	D0mass_Histo->SetTitle("Invariant Mass of the D0 ; Mass [GeV] ; Events ");
	D0mass_Histo->SetName("D0mass_Histo");
	if (debug) cout << "Teste 3 --------------------"<< endl;
	//End Histograms
	//-----------------------------------------------------------------------------

	//-------Reading the root file and the tree-------------------------------------	
	TFile *f1 = new TFile("D0DstarMC_25.root");
	TTree *t1 = (TTree*)f1->Get("analysis/data");
	if (debug)cout << "Teste 4 --------------------"<< endl;
	//---------------------------------------------------------------------------------
	// addressing the memory to vector and variables for file
		
	//For Vectors	
	TBranch *b_D0mass = t1->GetBranch("D0mass");
	b_D0mass->SetAddress(&D0mass);
	

	if (debug)cout << "Teste 5 --------------------"<< endl;

	//**********************************************************		
	//Reading Number of tree entries for file f1
	Long64_t nentries = t1->GetEntries();
	cout<< "Number of tree entries: "<< nentries <<std::endl;
	if (debug)cout << "Teste 6 --------------------"<< endl;
	//Long64_t GetEntriesFast = t1->GetEntriesFast();
	//cout<< "GetEntriesFast: "<< GetEntriesFast <<std::endl;

	for (Long64_t jentry=0; jentry < 100000; jentry++) //loop tree entries for file f1
	{
		if (debug)cout << "Teste 7 --------------------"<< endl;
		Long64_t ientry = t1->LoadTree(jentry);
      //std::cout << "nentries " << nentries << " ientry " << ientry << " jentry " << jentry <<std::endl;
      if (ientry < 0) break;
		
		if (debug)cout << "Teste 8 --------------------"<< endl;
		double percent = (jentry*100)/100000;
		if (jentry % 10000 == 0) cout <<"*****"<< percent << "per cent done***" << endl;
		//cout <<"*****entrada"<< jentry << "***" << endl;
		if (debug)cout << "Teste 9 --------------------"<< endl;
		
		b_D0mass->GetEntry(ientry);
		
		if (debug)cout << "Teste 10 --------------------"<< endl;
		//For D0
		for(unsigned i=0; i < D0mass->size(); i++)
		{  
			if (debug)cout << "Teste 11 --------------------"<< endl;
			D0mass_Histo->Fill(D0mass->at(i));
			//vectorInvariantMass_D0.push_back(D0mass->at(i));
			if (debug)cout << "Teste 12 --------------------"<< endl;
  		}
		
				if (debug)cout << "Teste 13 --------------------"<< endl;	
	}//End loop tree entries for file f1


	if (debug)cout << "Teste 14 --------------------"<< endl;
	D0mass_Histo->SetMarkerStyle(8);	
	D0mass_Histo->Sumw2();
	D0mass_Histo->Draw("e1p");

	if (debug)cout << "Teste 15 --------------------"<< endl;

	TFile f_analysis("DMesonsMassRegion.root","recreate"); //Creates root file
	TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	//t_analysis.Branch("vectorInvariantMass_Dstar",&vectorInvariantMass_Dstar); //Creates a branch
	//t_analysis.Branch("vectorInvariantMass_D0",&vectorInvariantMass_D0); //Creates a branch
	t_analysis.Fill();
	D0mass_Histo->Write();//Write() if a file is open, this function writes a root objectics on it.
	t_analysis.Write();  //Write in the root file

	cout << "-------END PROGRAM-------------"<< endl;
	
}//end program
