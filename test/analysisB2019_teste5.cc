// Draw histograms in a canvas straight from branches of two ntuples

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
#endif
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
void Drawhisto( string BRANCH, string NAME, string path)
{	//Creating Canvas
	
	TChain chain("analysis/data");
   chain.Add("D0DstarDataBparking_50.root");
   chain.Add("D0DstarDataBparking_54.root");
	//TH1F Histo;
	TCanvas* canvas = new TCanvas("canvas","",1200,600);
	//TH1::AddDirectory(kTRUE);
	chain.Draw(BRANCH.c_str()">>Histo");
	//chain.Print();
	//cout << "test21" << endl;
	//gDirectory->pwd();
	//gDirectory->ls();
	//cout << "test22" << endl;
	//TH1 *Histo = (TH1*)gDirectory->Get("Histo"); 
	//cout << NAME.c_str() << endl;
	//Histo->SetTitle(NAME.c_str());
	//canvas->Update();
	cout << Histo->GetTitle() << endl;
	//Histo->SetMarkerStyle(7);
	Histo->Draw("HIST");
	canvas->SaveAs(Form("%s%s.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void analysisB2019_teste5(){

	//call a file for a histogram style (Optional)
	//gROOT->LoadMacro("styleTDR.C"); 
	//setTDRStyle();

	string path = "/eos/user/r/ragomesd/analysisB2019/canvas/";
	
	//string save = Form("%s%s.png",path.c_str(),"15")
	string Branch = "D0pt";
	Drawhisto(Branch.c_str(),"Chi2 distribuition of the Kaon ; #Chi^{2} ; Events ", path.c_str());
	
}//end program
