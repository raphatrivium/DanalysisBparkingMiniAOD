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

int analysisB2019_teste()
{
	//call a file for a histogram style (Optional)
	gROOT->LoadMacro("styleTDR.C"); 
	setTDRStyle();
	
	//-------Reading the root file and the tree-------------------------------------	
	TFile *f1 = new TFile("D0DstarData_27.root");
	TTree *t1 = (TTree*)f1->Get("analysis/data");

	//**********************************************************		
	//Reading Number of tree entries for file f1
	Long64_t nentries = t1->GetEntries();
	cout<< "Number of tree entries: "<< nentries <<std::endl;

	//Long64_t GetEntriesFast = t1->GetEntriesFast();
	//cout<< "GetEntriesFast: "<< GetEntriesFast <<std::endl;

	/*TFile f_analysis("DMesonsMass.root","recreate"); //Creates root file
	TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	t_analysis.Branch("vectorInvariantMass_Dstar",&vectorInvariantMass_Dstar); //Creates a branch
	t_analysis.Branch("vectorInvariantMass_D0",&vectorInvariantMass_D0); //Creates a branch
	t_analysis.Fill();
	Dsmass_Histo->Write(); //Write() if a file is open, this function writes a root objectics on it.
	D0mass_Histo->Write(); //Write() if a file is open, this function writes a root objectics on it.
	t_analysis.Write();  //Write in the root file*/

    //========Creating Canvas=================================================================
	//========Creating Histgrams==============================================================
	TCanvas* c0 = new TCanvas("c0","Canvas 0 - behavior of the D",1200,600);
	c0->Divide(2);
	c0->cd(1);
	t1->Draw("D0mass>>D0massHisto(50,1.75,2.0)");
	D0massHisto->SetMarkerStyle(8);
	D0massHisto->Sumw2();
	D0massHisto->Draw("e1p");
	//a->GetPrimitive("h");
	D0massHisto->SetTitle("Invariant Mass of the D0 ; Mass [GeV] ; Events ");
	//--------------------------------------------------------------------------
	c0->cd(2);
	t1->Draw("Dsmass>>DsmassHisto(50,1.85,2.15)");
	DsmassHisto->SetMarkerStyle(8);
	DsmassHisto->Sumw2();
	DsmassHisto->Draw("e1p");
	//a->GetPrimitive("h");
	DsmassHisto->SetTitle("Invariant Mass of the D* ; Mass [GeV] ; Events ");
	c0->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas0.png");

	//=========================================================================
	
	TCanvas* c1 = new TCanvas("c1","Canvas 1 - behavior of the D",1200,600);
	c1->Divide(2);
	c1->cd(1);
	t1->Draw("D0pt>>D0ptHisto(100,0,40)");
	D0ptHisto->Draw("histo");
	TLegend* leg_D0mass = new TLegend(0.5,0.5,0.75,0.65);
   	leg_D0mass->SetFillColor(kWhite);
	leg_D0mass->SetFillStyle(1001);
	leg_D0mass->SetBorderSize(0);
	leg_D0mass->AddEntry(D0ptHisto,"pT > 3 GeV","SAME");
	leg_D0mass->AddEntry(D0ptHisto,"Bparking1-6","SAME");
	leg_D0mass->Draw();
	//a->GetPrimitive("h");
	D0ptHisto->SetTitle("pT distribuition of the D0 ; pT [GeV] ; Events ");
	//--------------------------------------------------------------------------
	c1->cd(2);
	t1->Draw("Dspt>>DsptHisto(100,0,40)");
	DsptHisto->Draw("histo");
	//a->GetPrimitive("h");
	DsptHisto->SetTitle("pT distribuition of the D* ; pT [GeV] ; Events ");
	c1->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas1.png");
	//=========================================================================
	
	TCanvas* c2 = new TCanvas("c2","Canvas 2 - behavior of the D",1200,600);
	c2->Divide(2);
	c2->cd(1);
	t1->Draw("D0eta>>D0etaHisto(100,-3,3)");
	D0etaHisto->Draw("histo");
	//a->GetPrimitive("h");
	D0etaHisto->SetTitle("Pseudo-rapidity distribuition of the D0 ; #eta ; Events ");
	//--------------------------------------------------------------------------
	c2->cd(2);
	t1->Draw("Dseta>>DsetaHisto(100,-3,3)");
	DsetaHisto->Draw("histo");
	//a->GetPrimitive("h");
	DsetaHisto->SetTitle("Pseudo-rapidity distribuition of the D* ; #eta ; Events ");
	c2->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas2.png");
	//=========================================================================

	TCanvas* c31 = new TCanvas("c31","Canvas c31 - behavior of the D",1200,600);
    t1->Draw("fabs(D0eta):D0pt>>D0etaD0ptHisto(100,0,40,100,0,4)");
   	D0etaD0ptHisto->Draw("CONTZ");
	//D0etaD0ptHisto->Draw("CONT4COLZ")
   	//a->GetPrimitive("h");
   	D0etaD0ptHisto->SetTitle(" Pseudo-rapidity x Tranverse Momentum D0 ; pT [GeV] ; #Eta ");
    c31->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas31.png");
	//=========================================================================
	
	TCanvas* c3 = new TCanvas("c3","Canvas 3 - behavior of the D",1200,600);
	c3->Divide(2);
	c3->cd(1);
	t1->Draw("D0phi>>D0phiHisto(100,-4,4)");
	D0phiHisto->Draw("histo");
	//a->GetPrimitive("h");
	D0phiHisto->SetTitle("Phi Angle distribuition of the D0 ; #phi ; Events ");
	//--------------------------------------------------------------------------
	c3->cd(2);
	t1->Draw("Dsphi>>DsphiHisto(100,-4,4)");
	DsphiHisto->Draw("histo");
	//a->GetPrimitive("h");
	DsphiHisto->SetTitle("Phi Angle distribuition of the D* ; #phi ; Events ");
	c3->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas3.png");
	//=========================================================================
	
	TCanvas* c4 = new TCanvas("c4","Canvas 4 - behavior of the SlowPion",1200,600);
	c4->Divide(2);
	c4->cd(1);
	t1->Draw("TrkSpt>>TrkSptHisto(100,0,2.5)");
	TrkSptHisto->Draw("histo");
	TLegend* leg_TrkSpt = new TLegend(0.5,0.5,0.75,0.65);
   	leg_TrkSpt->SetFillColor(kWhite);
	leg_TrkSpt->SetFillStyle(1001);
	leg_TrkSpt->SetBorderSize(0);
	leg_TrkSpt->AddEntry(TrkSptHisto,"pT > 0.3 GeV","SAME");
	leg_TrkSpt->AddEntry(TrkSptHisto,"Bparking1-6","SAME");
	leg_TrkSpt->Draw();
	//a->GetPrimitive("h");
	TrkSptHisto->SetTitle("pt distribuition of the SlowPion ; pT [GeV] ; Events ");
	//--------------------------------------------------------------------------
		
	c4->cd(2);
	t1->Draw("TrkSnhits>>TrkSnhitsHisto(100,0,45)");
	TrkSnhitsHisto->Draw("histo");
	//a->GetPrimitive("h");
	TLegend* leg_TrkSnhits = new TLegend(0.6,0.5,0.85,0.65);
   	leg_TrkSnhits->SetFillColor(kWhite);
	leg_TrkSnhits->SetFillStyle(1001);
	leg_TrkSnhits->SetBorderSize(0);
	leg_TrkSnhits->AddEntry(TrkSnhitsHisto,"Nhits > 2","SAME");
	leg_TrkSnhits->AddEntry(TrkSnhitsHisto,"Bparking1-6","SAME");
	leg_TrkSnhits->Draw();
	TrkSnhitsHisto->SetTitle("nhits distribuition of the SlowPion ; Number of hits ; Events ");
	c4->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas4.png");
	//=========================================================================

	TCanvas* c5 = new TCanvas("c5","Canvas 5 - behavior of the SlowPion",1200,600);
	c5->Divide(2);
	c5->cd(1);
	t1->Draw("TrkSdxy>>TrkSdxyHisto(100,-1.5,1.5)");
	TrkSdxyHisto->Draw("histo");
	//a->GetPrimitive("h");
	TLegend* leg_TrkSdxy = new TLegend(0.2,0.5,0.45,0.65);
   	leg_TrkSdxy->SetFillColor(kWhite);
	leg_TrkSdxy->SetFillStyle(1001);
	leg_TrkSdxy->SetBorderSize(0);
	leg_TrkSdxy->AddEntry(TrkSdxyHisto,"Dxy<3cm","SAME");
	leg_TrkSdxy->AddEntry(TrkSdxyHisto,"Bparking1-6","SAME");
	leg_TrkSdxy->Draw();
	TrkSdxyHisto->SetTitle("dxy distribuition of the SlowPion ; dx [cm] ; Events ");
	//--------------------------------------------------------------------------
	c5->cd(2);
	t1->Draw("TrkSdz>>TrkSdzHisto(100,-4,4)");
	TrkSdzHisto->Draw("histo");
	//a->GetPrimitive("h");
	TLegend* leg_TrkSdz = new TLegend(0.2,0.5,0.45,0.65);
   	leg_TrkSdz->SetFillColor(kWhite);
	leg_TrkSdz->SetFillStyle(1001);
	leg_TrkSdz->SetBorderSize(0);
	leg_TrkSdz->AddEntry(TrkSdzHisto,"Dz<3cm","SAME");
	leg_TrkSdz->AddEntry(TrkSdzHisto,"Bparking1-6","SAME");
	leg_TrkSdz->Draw();
	TrkSdzHisto->SetTitle("dz distribuition of the SlowPion ; dz [cm] ; Events ");
	c5->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas5.png");
	//=========================================================================

	TCanvas* c6 = new TCanvas("c6","Canvas 6 - behavior of the SlowPion",1200,600);
	c6->Divide(2);
	c6->cd(1);
	t1->Draw("TrkSeta>>TrkSetaHisto(100,-3.1,3.1)");
	TrkSetaHisto->Draw("histo");
	//a->GetPrimitive("h");
	TLegend* leg_TrkSeta = new TLegend(0.4,0.2,0.65,0.45);
   	leg_TrkSeta->SetFillColor(kWhite);
	leg_TrkSeta->SetFillStyle(1001);
	leg_TrkSeta->SetBorderSize(0);
	leg_TrkSeta->AddEntry(TrkSetaHisto,"#eta<3 ","SAME");
	leg_TrkSeta->AddEntry(TrkSetaHisto,"Bparking1-6","SAME");
	leg_TrkSeta->Draw();
	TrkSetaHisto->SetTitle("eta distribuition of the SlowPion ; #eta ; Events ");
	//--------------------------------------------------------------------------
	c6->cd(2);
	t1->Draw("TrkSphi>>TrkSphiHisto(100,-4,4)");
	TrkSphiHisto->Draw("histo");
	//a->GetPrimitive("h");
	TrkSphiHisto->SetTitle("Phi distribuition of the SlowPion ; #Phi ; Events ");
	c6->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas6.png");
	//=========================================================================

	TCanvas* c7 = new TCanvas("c7","Canvas 7 - behavior of the SlowPion",1200,600);
	c7->Divide(2);
	c7->cd(1);
	t1->Draw("TrkSchi2>>TrkSchi2Histo(100,0,3)");
	TrkSchi2Histo->Draw("histo");
	//a->GetPrimitive("h");
	TrkSchi2Histo->SetTitle("Chi2 distribuition of the SlowPion ; #Chi^{2} ; Events ");
	//--------------------------------------------------------------------------
	//c6->cd(2);
	


	c7->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas7.png");
	//=========================================================================

	TCanvas* c8 = new TCanvas("c8","Canvas 8 - behavior of the Pion",1200,600);
	c8->Divide(2);
	c8->cd(1);
	t1->Draw("Trkpipt>>TrkpiptHisto(100,0,3.)");
	TrkpiptHisto->Draw("histo");
	TLegend* leg_Trkpipt = new TLegend(0.5,0.5,0.75,0.65);
	leg_Trkpipt->SetFillColor(kWhite);
	leg_Trkpipt->SetFillStyle(1001);
	leg_Trkpipt->SetBorderSize(0);
	leg_Trkpipt->AddEntry(TrkpiptHisto,"pT>0.5 GeV","SAME");
	leg_Trkpipt->AddEntry(TrkpiptHisto,"Bparking1-6","SAME");
	leg_Trkpipt->Draw();
	//a->GetPrimitive("h");
	TrkpiptHisto->SetTitle("pt distribuition of the Pion ; pT [GeV] ; Events ");
	//--------------------------------------------------------------------------
	c8->cd(2);
	t1->Draw("Trkpinhits>>TrkpinhitsHisto(100,0,45)");
	TrkpinhitsHisto->Draw("histo");
	//a->GetPrimitive("h");
	TLegend* leg_Trkpinhits = new TLegend(0.6,0.5,0.85,0.65);
   	leg_Trkpinhits->SetFillColor(kWhite);
	leg_Trkpinhits->SetFillStyle(1001);
	leg_Trkpinhits->SetBorderSize(0);
	leg_Trkpinhits->AddEntry(TrkpinhitsHisto,"Nhits #geq 5","SAME");
	leg_Trkpinhits->AddEntry(TrkpinhitsHisto,"Bparking1-6","SAME");
	leg_Trkpinhits->Draw();
	TrkpinhitsHisto->SetTitle("nhits distribuition of the Pion ; Number of hits ; Events ");
	c8->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas8.png");
	//=========================================================================

	TCanvas* c9 = new TCanvas("c9","Canvas 9 - behavior of the Pion",1200,600);
	c9->Divide(2);
	c9->cd(1);
	t1->Draw("Trkpidxy>>TrkpidxyHisto(100,-1,1)");
	TrkpidxyHisto->Draw("histo");
	//a->GetPrimitive("h");
	TLegend* leg_Trkpidxy = new TLegend(0.2,0.5,0.45,0.65);
   	leg_Trkpidxy->SetFillColor(kWhite);
	leg_Trkpidxy->SetFillStyle(1001);
	leg_Trkpidxy->SetBorderSize(0);
	leg_Trkpidxy->AddEntry(TrkpidxyHisto,"Dxy<0.1 cm","SAME");
	leg_Trkpidxy->AddEntry(TrkpidxyHisto,"Bparking1-6","SAME");
	leg_Trkpidxy->Draw();
	TrkpidxyHisto->SetTitle("dxy distribuition of the Pion ; dx [cm] ; Events ");
	//--------------------------------------------------------------------------
	c9->cd(2);
	t1->Draw("Trkpidz>>TrkpidzHisto(100,-1.5,1.5)");
	TrkpidzHisto->Draw("histo");
	//a->GetPrimitive("h");
	TLegend* leg_Trkpidz = new TLegend(0.2,0.5,0.45,0.65);
   	leg_Trkpidz->SetFillColor(kWhite);
	leg_Trkpidz->SetFillStyle(1001);
	leg_Trkpidz->SetBorderSize(0);
	leg_Trkpidz->AddEntry(TrkpidzHisto,"Dz<1 cm","SAME");
	leg_Trkpidz->AddEntry(TrkpidzHisto,"Bparking1-6","SAME");
	leg_Trkpidz->Draw();
	TrkpidzHisto->SetTitle("dz distribuition of the Pion ; dz [cm] ; Events ");
	c9->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas9.png");
	//=========================================================================

	TCanvas* c10 = new TCanvas("c10","Canvas 10 - behavior of the Pion",1200,600);
	c10->Divide(2);
	c10->cd(1);
	t1->Draw("Trkpieta>>TrkpietaHisto(100,-3,3)");
	TrkpietaHisto->Draw("histo");
	//a->GetPrimitive("h");
	TLegend* leg_Trkpieta = new TLegend(0.4,0.2,0.65,0.45);
   	leg_Trkpieta->SetFillColor(kWhite);
	leg_Trkpieta->SetFillStyle(1001);
	leg_Trkpieta->SetBorderSize(0);
	leg_Trkpieta->AddEntry(TrkpietaHisto,"#eta<3 ","SAME");
	leg_Trkpieta->AddEntry(TrkpietaHisto,"Bparking1-6","SAME");
	leg_Trkpieta->Draw();
	TrkpietaHisto->SetTitle("eta distribuition of the Pion ; #eta ; Events ");
	//--------------------------------------------------------------------------
	c10->cd(2);
	t1->Draw("Trkpiphi>>TrkpiphiHisto(100,-4,4)");
	TrkpiphiHisto->Draw("HIST");
	//a->GetPrimitive("h");
	TrkpiphiHisto->SetTitle("Phi distribuition of the Pion ; #Phi ; Events ");
	c10->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas10.png");
	//=========================================================================

	TCanvas* c11 = new TCanvas("c11","Canvas 11 - behavior of the Pion",1200,600);
	c11->Divide(2);
	c11->cd(1);
	t1->Draw("Trkpichi2>>Trkpichi2Histo(100,0,3)");
	Trkpichi2Histo->Draw("histo");
	//a->GetPrimitive("h");
	Trkpichi2Histo->SetTitle("Chi2 distribuition of the Pion ; #Chi^{2} ; Events ");
	//--------------------------------------------------------------------------
	//c11->cd(2);
	


	c11->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas11.png");
	//=========================================================================

	TCanvas* c12 = new TCanvas("c12","Canvas 12 - behavior of the Kaon",1200,600);
	c12->Divide(2);
	c12->cd(1);
	t1->Draw("TrkKpt>>TrkKptHisto(100,0,3.)");
	TrkKptHisto->Draw("HIST");
	//a->GetPrimitive("h");
	TLegend* leg_TrkKpt = new TLegend(0.5,0.5,0.75,0.65);
	leg_TrkKpt->SetFillColor(kWhite);
	leg_TrkKpt->SetFillStyle(1001);
	leg_TrkKpt->SetBorderSize(0);
	leg_TrkKpt->AddEntry(TrkKptHisto,"pT>0.5 GeV","SAME");
	leg_TrkKpt->AddEntry(TrkKptHisto,"Bparking1-6","SAME");
	leg_TrkKpt->Draw();
	TrkKptHisto->SetTitle("pt distribuition of the Kaon ; pT [GeV] ; Events ");
	//--------------------------------------------------------------------------
	c12->cd(2);
	t1->Draw("TrkKnhits>>TrkKnhitsHisto(100,0,45)");
	TrkKnhitsHisto->Draw("HIST");
	//a->GetPrimitive("h");
	TLegend* leg_TrkKnhits = new TLegend(0.6,0.5,0.85,0.65);
   	leg_TrkKnhits->SetFillColor(kWhite);
	leg_TrkKnhits->SetFillStyle(1001);
	leg_TrkKnhits->SetBorderSize(0);
	leg_TrkKnhits->AddEntry(TrkKnhitsHisto,"Nhits #geq 5","SAME");
	leg_TrkKnhits->AddEntry(TrkKnhitsHisto,"Bparking1-6","SAME");
	leg_TrkKnhits->Draw();
	TrkKnhitsHisto->SetTitle("nhits distribuition of the Kaon ; Number of hits ; Events ");
	c12->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas12.png");
	//=========================================================================

	TCanvas* c13 = new TCanvas("c13","Canvas 12 - behavior of the Kaon",1200,600);
	c13->Divide(2);
	c13->cd(1);
	t1->Draw("TrkKdxy>>TrkKdxyHisto(100,-1,1)");
	TrkKdxyHisto->Draw("HIST");
	//a->GetPrimitive("h");
	TLegend* leg_TrkKdxy = new TLegend(0.2,0.5,0.45,0.65);
   	leg_TrkKdxy->SetFillColor(kWhite);
	leg_TrkKdxy->SetFillStyle(1001);
	leg_TrkKdxy->SetBorderSize(0);
	leg_TrkKdxy->AddEntry(TrkKdxyHisto,"Dxy<0.1 cm","SAME");
	leg_TrkKdxy->AddEntry(TrkKdxyHisto,"Bparking1-6","SAME");
	leg_TrkKdxy->Draw();
	TrkKdxyHisto->SetTitle("dxy distribuition of the Kaon ; dx [cm] ; Events ");
	//--------------------------------------------------------------------------
	c13->cd(2);
	t1->Draw("TrkKdz>>TrkKdzHisto(100,-1.5,1.5)");
	TrkKdzHisto->Draw("HIST");
	//a->GetPrimitive("h");
	TLegend* leg_TrkKdz = new TLegend(0.2,0.5,0.45,0.65);
   	leg_TrkKdz->SetFillColor(kWhite);
	leg_TrkKdz->SetFillStyle(1001);
	leg_TrkKdz->SetBorderSize(0);
	leg_TrkKdz->AddEntry(TrkKdzHisto,"Dz<1 cm","SAME");
	leg_TrkKdz->AddEntry(TrkKdzHisto,"Bparking1-6","SAME");
	leg_TrkKdz->Draw();
	TrkKdzHisto->SetTitle("dz distribuition of the Pion ; dz [cm] ; Events ");
	c13->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas13.png");
	//=========================================================================

	TCanvas* c14 = new TCanvas("c14","Canvas 14 - behavior of the Kaon",1200,600);
	c14->Divide(2);
	c14->cd(1);
	t1->Draw("TrkKeta>>TrkKetaHisto(100,-3,3)");
	TrkKetaHisto->Draw("HIST");
	//a->GetPrimitive("h");
	TLegend* leg_TrkKeta = new TLegend(0.4,0.2,0.65,0.45);
   	leg_TrkKeta->SetFillColor(kWhite);
	leg_TrkKeta->SetFillStyle(1001);
	leg_TrkKeta->SetBorderSize(0);
	leg_TrkKeta->AddEntry(TrkKetaHisto,"#eta<3 ","SAME");
	leg_TrkKeta->AddEntry(TrkKetaHisto,"Bparking1-6","SAME");
	leg_TrkKeta->Draw();
	TrkKetaHisto->SetTitle("eta distribuition of the Kaon ; #eta ; Events ");
	//--------------------------------------------------------------------------
	c14->cd(2);
	t1->Draw("TrkKphi>>TrkKphiHisto(100,-4,4)");
	TrkKphiHisto->Draw("HIST");
	//a->GetPrimitive("h");
	TrkKphiHisto->SetTitle("Phi distribuition of the Kaon ; #Phi ; Events ");
	c14->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas14.png");
	//=========================================================================

	TCanvas* c15 = new TCanvas("c15","Canvas 15 - behavior of the Kaon",1200,600);
	c15->Divide(2);
	c15->cd(1);
	t1->Draw("TrkKchi2>>TrkKchi2Histo(100,0,3)");
	TrkKchi2Histo->Draw("HIST");
	//a->GetPrimitive("h");
	TrkKchi2Histo->SetTitle("Chi2 distribuition of the Kaon ; #Chi^{2} ; Events ");
	//--------------------------------------------------------------------------
	//c6->cd(2);
	


	c15->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas15.png");
	//=========================================================================

	TCanvas* c16 = new TCanvas("c16","Canvas 16 - behavior of the Kaon",1200,600);
	//c16->Divide(2);
	//c16->cd(1);
	t1->Draw("D0_VtxProb>>D0VtxProbHisto(100,0,0.1)");
	D0VtxProbHisto->Draw("HIST");
	//a->GetPrimitive("h");
	D0VtxProbHisto->SetTitle("SV Confidence Level ; CL ; Events ");
	c16->SaveAs("/eos/user/r/ragomesd/analysisB2019/canvas/Canvas16.png");
	//--------------------------------------------------------------------------



	/*TFile f_analysis("Ds_D0_K_pi_Bparking.root","recreate"); //Creates root file
	TTree t_analysis("t_analysis","analise_Tree"); //Creates a Tree
	//t_analysis.Branch("vectorInvariantMass_Dstar",&vectorInvariantMass_Dstar); //Creates a branch
	//t_analysis.Branch("vectorInvariantMass_D0",&vectorInvariantMass_D0); //Creates a branch
	//t_analysis.Fill();
	D0massHisto.Write(); //Write() if a file is open, this function writes a root objectics on it.
	D0ptHisto->Write();
	D0etaHisto->Write();
	D0phiHisto->Write();
	DsmassHisto->Write();
	DsptHisto->Write();
	DsetaHisto->Write();
	DsphiHisto->Write();
	D0etaD0ptHisto->Write();
	TrkSptHisto->Write();
	TrkSnhitsHisto->Write();
	TrkSdxyHisto->Write();
	TrkSdzHisto->Write();
	TrkSphiHisto->Write();
	TrkSetaHisto->Write();
	TrkSchi2Histo->Write();
	TrkpiptHisto->Write();
	TrkpinhitsHisto->Write();
	TrkpidxyHisto->Write();
	TrkpidzHisto->Write();
	TrkpietaHisto->Write();
	TrkpiphiHisto->Write();
	Trkpichi2Histo->Write();
	TrkKptHisto->Write();
	TrkKnhitsHisto->Write();
	TrkKdxyHisto->Write();
	TrkKdzHisto->Write();
	TrkKetaHisto->Write();
	TrkKphiHisto->Write();
	TrkKchi2Histo->Write();
	t_analysis.Write();  //Write in the root file*/
	

	
	



		
}//end program

