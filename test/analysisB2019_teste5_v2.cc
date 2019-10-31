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
void Drawhisto( string BRANCH, string NBIN , string min , string max, string NAME, string path)
{	//Creating Canvas
	
	TChain chain("analysis/data");
   //chain.Add("D0DstarDataBparking_50.root");
   //chain.Add("D0DstarDataBparking_54.root");
	chain.Add("merge_data1.root");
	chain.Add("merge_data2.root");
	chain.Add("merge_data3.root");
	chain.Add("merge_data4.root");
	chain.Add("merge_data5.root");
	chain.Add("merge_data6.root");
	
	TCanvas* canvas = new TCanvas("canvas","",1200,600);
	chain.Draw(Form("%s>>Histo(%s,%s,%s)",BRANCH.c_str(), NBIN.c_str(), min.c_str(), max.c_str()));
	TH1 *Teste = (TH1*)gDirectory->Get("Histo");
	Teste->SetTitle(NAME.c_str());
	Teste->SetMarkerStyle(7);
	Teste->SetFillColor(kRed);
	Teste->Draw("HIST");
	
	canvas->SaveAs(Form("%s%s.png", path.c_str(), BRANCH.c_str()));
	delete canvas;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void analysisB2019_teste5(){

	//call a file for a histogram style (Optional)
	gROOT->LoadMacro("styleTDR.C"); 
	setTDRStyle();

	//Read varous files as one (must be in same structure and names)
	TChain chain("analysis/data");
   //chain.Add("D0DstarDataBparking_50.root");
   //chain.Add("D0DstarDataBparking_54.root");
	chain.Add("merge_data1.root");
	chain.Add("merge_data2.root");
	chain.Add("merge_data3.root");
	chain.Add("merge_data4.root");
	chain.Add("merge_data5.root");
	chain.Add("merge_data6.root");

	//Creates root file and Tree
	TFile f_analysis("Ds_D0_K_pi_Bparking.root","recreate"); 
	TTree t_analysis("t_analysis","analise_Tree");

	//========Creating Canvas=================================================================
	TCanvas* c900 = new TCanvas("c900","Canvas c900 - behavior of the D",1200,600);
	chain.Draw("fabs(D0eta):D0pt>>D0etaD0ptHisto(100,3,15,100,0,1.5)");
	D0etaD0ptHisto->Draw("COLZ"); //COL,COLZ,CONT0,CONTZ,CONT4COLZ
   //a->GetPrimitive("h");
   D0etaD0ptHisto->SetTitle("Pseudo-rapidity x Tranverse Momentum D0 ; pT [GeV] ; #eta ");
	c900->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/TH2D0etaPT.png");
	
	//========Creating Canvas==============================================================
	TCanvas* c0 = new TCanvas("c0","Canvas 0 - behavior of the D",1200,600);
	c0->Divide(2);
	c0->cd(1);
	chain.Draw("D0mass>>D0massHistoTemp");
	TH1 *D0massHisto = (TH1*)gDirectory->Get("D0massHistoTemp");
	D0massHisto->SetTitle("Invariant Mass of the D0 ; Mass [GeV] ; Events ");
	D0massHisto->SetMarkerStyle(7);
	D0massHisto->Sumw2();
	D0massHisto->Draw("e1p");

	chain.Draw("Dsmass>>DsmassHistoTemp");
	TH1 *DsmassHisto = (TH1*)gDirectory->Get("DsmassHistoTemp");
	DsmassHisto->SetTitle("Invariant Mass of the D* ; Mass [GeV] ; Events ");
	c0->cd(2);
	DsmassHisto->SetLineColor(kBlue);
	DsmassHisto->SetMarkerStyle(7);
	DsmassHisto->Sumw2();
	DsmassHisto->Draw("e1p");
	
	c0->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/MassD0Ds.png");

	//---------------------------------------------------------------------------------------------

	//local to save the png's
	string path = "/eos/user/r/ragomesd/analysisB2019/canvas/";
	
	//string save = Form("%s%s.png",path.c_str(),"15")
	string Branch = "D0pt"; string nbin = "100"; string nmin = "0"; string nmax = "30";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"pT distribuition of the D0 ; pT [GeV] ; Events ", path.c_str());
	string Branch = "Dspt"; string nbin = "100"; string nmin = "0"; string nmax = "30";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"pT distribuition of the D* ; pT [GeV] ; Events ", path.c_str());
	string Branch = "D0eta"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"Pseudo-rapidity distribuition of the D0 ; #eta ; Events", path.c_str());
	string Branch = "Dseta"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"Pseudo-rapidity distribuition of the D* ; #eta ; Events ", path.c_str());
	string Branch = "D0phi"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"#Phi distribuition of the D0 ; #Phi ; Events ", path.c_str());
	string Branch = "Dsphi"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"#Phi distribuition of the D* ; #Phi ; Events ", path.c_str());
	string Branch = "TrkSpt"; string nbin = "100"; string nmin = "0"; string nmax = "30";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"pt distribuition of the SlowPion ; pT [GeV] ; Events ", path.c_str());	
	string Branch = "TrkSnhits"; string nbin = "100"; string nmin = "0"; string nmax = "50";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"nhits distribuition of the SlowPion ; Number of hits ; Events ", path.c_str());
	string Branch = "TrkSdxy"; string nbin = "100"; string nmin = "-0.4"; string nmax = "0.4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dxy distribuition of the SlowPion ; dx [cm] ; Events ", path.c_str());
	string Branch = "TrkSdz"; string nbin = "100"; string nmin = "-0.4"; string nmax = "0.4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dz distribuition of the SlowPion ; dz [cm] ; Events ", path.c_str());
	string Branch = "TrkSeta"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"eta distribuition of the SlowPion ; #eta ; Events ", path.c_str());
	string Branch = "TrkSphi"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"Phi distribuition of the SlowPion ; #Phi ; Events ", path.c_str());
	string Branch = "TrkSchi2"; string nbin = "100"; string nmin = "0"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"Chi2 distribuition of the SlowPion ; #Chi^{2} ; Events ", path.c_str());
	string Branch = "Trkpipt"; string nbin = "100"; string nmin = "0"; string nmax = "30";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"pt distribuition of the Pion ; pT [GeV] ; Events ", path.c_str());
	string Branch = "Trkpinhits"; string nbin = "100"; string nmin = "0"; string nmax = "50";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"nhits distribuition of the Pion ; Number of hits ; Events ", path.c_str());
	string Branch = "Trkpidxy"; string nbin = "100"; string nmin = "-0.4"; string nmax = "0.4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dxy distribuition of the Pion ; dx [cm] ; Events ", path.c_str());
	string Branch = "Trkpidz"; string nbin = "100"; string nmin = "-0.4"; string nmax = "0.4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dz distribuition of the Pion ; dz [cm] ; Events ", path.c_str());
	string Branch = "Trkpieta"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"eta distribuition of the Pion ; #eta ; Events ", path.c_str());
	string Branch = "Trkpiphi"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"phi distribuition of the Pion ; #Phi ; Events ", path.c_str());
	string Branch = "Trkpichi2"; string nbin = "100"; string nmin = "0"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"Chi2 distribuition of the Pion ; #Chi^{2} ; Events ", path.c_str());
	string Branch = "TrkKpt"; string nbin = "100"; string nmin = "0"; string nmax = "30";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"pt distribuition of the Kaon ; pT [GeV] ; Events ", path.c_str());
	string Branch = "TrkKnhits"; string nbin = "100"; string nmin = "0"; string nmax = "50";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"nhits distribuition of the Kaon ; Number of hits ; Events ", path.c_str());
	string Branch = "TrkKdxy"; string nbin = "100"; string nmin = "-0.4"; string nmax = "0.4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dxy distribuition of the Kaon ; dx [cm] ; Events ", path.c_str());
	string Branch = "TrkKdz"; string nbin = "100"; string nmin = "-0.4"; string nmax = "0.4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"dz distribuition of the Kaon ; dz [cm] ; Events ", path.c_str());
	string Branch = "TrkKeta"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"eta distribuition of the Kaon ; #eta ; Events ", path.c_str());
	string Branch = "TrkKphi"; string nbin = "100"; string nmin = "-4"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"phi distribuition of the Kaon ; #Phi ; Events ", path.c_str());
	string Branch = "TrkKchi2"; string nbin = "100"; string nmin = "0"; string nmax = "4";
	Drawhisto(Branch.c_str(), nbin.c_str(), nmin.c_str(), nmax.c_str(),"Chi2 distribuition of the Kaon ; #Chi^{2} ; Events ", path.c_str());

	D0massHisto->Write();
	DsmassHisto->Write();
	t_analysis.Write();  //Write in the root file*/
	//f_analysis.Close();
	
}//end program
