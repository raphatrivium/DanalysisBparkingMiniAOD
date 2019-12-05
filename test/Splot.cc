//PDF - Probability density function
#include <iomanip>
#include <sstream>
#include <vector>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TSystem.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooBifurGauss.h>
#include "TMath.h"
#include <RooGenericPdf.h>
#include "TRatioPlot.h"
#include <RooBifurGauss.h>
#include <RooProduct.h>
#include <RooHist.h>
#include "RooStats/SPlot.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include <iostream>
#include <TF1.h>
#include <RooPolynomial.h>
#include <fstream>
#include <TGraph.h>
#include "TMultiGraph.h"
using namespace RooStats;
using namespace RooFit;
using namespace std;

// see below for implementation
void set_up_workspace_variables(RooWorkspace& ws);
void AddData(RooWorkspace& ws);
void AddModel(RooWorkspace& ws);
void plot_complete_fit(RooWorkspace& ws);
void DoSPlot(RooWorkspace& ws);
void MakePlots(RooWorkspace& ws);

void DstarD0Fit_RooDataSet_v4Workspace()
{
	// Create a new workspace to manage the project.
  	RooWorkspace* MyWS = new RooWorkspace("MyWS");
	
	cout << "set up variables" << endl;
  	set_up_workspace_variables(*MyWS);

	// add some toy data to the workspace
	cout << "A D D   D A T A" << endl;	
	AddData(*MyWS);

	// Inside this function you will find a discription our model.
	cout << "A D D   M O D E L" << endl;
	AddModel(*MyWS);

	//Fit Signal + Background.
	plot_complete_fit(*MyWS);

	// do sPlot.This wil make a new dataset with sWeights added for every event.
	DoSPlot(*MyWS);

	// Make some plots showing the discriminating variable and
	// the control variable after unfolding.
	MakePlots(*MyWS);

	// inspect the workspace if you wish
	//MyWS->Print();

  	// cleanup
  	delete MyWS;
	
	cout << "---------------------------------" << endl;
	cout << "E N D   P R O G R A M" << endl;
	cout << "---------------------------------" << endl;

}//end program
//=====================================================================
//=====================================================================
//=====================================================================
void set_up_workspace_variables(RooWorkspace& ws)
{
	// set range of observable
	Double_t lowRange = 1.77, highRange = 1.95;
	//Double_t lowRange = 1.93, highRange = 2.10;

	// make a RooRealVar for the observables
	RooRealVar Dsmass("Dsmass","Dsmass M_{inv}",lowRange,highRange);
	//RooRealVar Dsmass("Dsmass","Dsmass",1.93,2.10);
	RooRealVar D0pt("D0pt","D0pt",0.,60.,"GeV");

	RooRealVar DsmassWrong("DsmassWrong","DsmassWrong",0.,30.,"GeV");

	ws.import(Dsmass);
	ws.import(D0pt);
	ws.import(DsmassWrong);

}
//=====================================================================
void AddData(RooWorkspace& ws){
	// Add a dataset
	TFile* file = new TFile("DMesonsHistograms_Data.root");
	TTree* t1_data = (TTree *)file->Get("t_analysis");
	cout << "Reading root file" << endl;

	RooArgList arg_list ("arg_list");

	arg_list.add(*(ws.var("Dsmass")));
	arg_list.add(*(ws.var("D0pt")));

	RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);

	ws.import(*data, Rename("data"));

}//End addData
//=====================================================================
void AddModel(RooWorkspace& ws) {

	RooRealVar Dsmass = *(ws.var("Dsmass"));
	RooDataSet* data = (RooDataSet*) ws.data("data");

	///------------------------------------------------
	// P D F   S I G N A L 
	// ------------------------------------------------
	//RooRealVar mean( "mean", "mean", 1.864, 1.840, 1.900);
	RooRealVar mean( "mean", "mean", 2.010, 1.980, 2.100);
	//RooRealVar sigma1( "sigma1", "sigma1", 0.01748, 0.00, 1.0);
	//RooRealVar sigma2( "sigma2", "sigma2", 0.00924, 0.001, 1.);
	RooRealVar sigma1( "sigma1", "sigma1", 0.01740, 0.00, 1.0);
	RooRealVar sigma2( "sigma2", "sigma2", 0.00912, 0.001, 1.);
	RooGaussian gauss1("gauss1","gauss1",Dsmass,mean,sigma1);
	RooGaussian gauss2("gauss2","gauss2",Dsmass,mean,sigma2);
	RooRealVar coef("coef","fraction of Gaussian1", 0.452, 0.0, 1.0);
	//RooRealVar coef("coef","fraction of Gaussian1", 1);
	RooAddPdf  sig("sig","gauss1 + gauss2",RooArgList(gauss1,gauss2),RooArgList(coef));

   sigma1.setConstant();
   sigma2.setConstant();
	
	cout << "Roofit variables ok" << endl;

	// Make the RooRealVar mean constant. It will not vary in the input range
	//mean.setConstant();
			
	//------------------------------------------------
	// P D F   B A C K G R O U N D
	// ------------------------------------------------
	//ROOFIT third degree polynomial
	RooRealVar lambda("lambda","lambda",-2.26,-10.,0.0);
  	RooExponential bkg("bkg", "bkg", Dsmass, lambda);
	
	//------------------------------------------------
	// A d d  s i g n a l   a n d   b a c k g r o u n d
	// ------------------------------------------------
	double n_signal_initial = data->sumEntries(TString::Format("abs(Dsmass-%g)<0.05",mean.getVal())) - data->sumEntries(TString::Format("abs(Dsmass-%g)<0.10&&abs(Dsmass-%g)>0.05",mean.getVal(),mean.getVal()));
  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
	// Sum the composite signal and background 
	//RooRealVar n_signal("n_signal","fraction of signal", 1200, 10, 2500);
	//RooRealVar n_combinatorial("n_combinatorial","fraction of background", 500, 110, 2500);
	RooAddPdf  model("model","sig + bkg", RooArgList(sig,bkg), RooArgList(n_signal,n_combinatorial));
	cout << "Roofit Parameters ok" << endl;

	ws.import(model);
}//End AddModel
//=====================================================================
void plot_complete_fit(RooWorkspace& ws){

	RooAbsPdf* model = ws.pdf("model");
	RooDataSet* data = (RooDataSet*) ws.data("data");

	RooRealVar Dsmass = *(ws.var("Dsmass"));
	RooRealVar* sigma1 = ws.var("sigma1");
	RooRealVar* sigma2 = ws.var("sigma2");

	RooRealVar DsmassWrong = *(ws.var("DsmassWrong"));
	

	model->fitTo(*data,Range("all"));

	//FIT
	RooPlot* massframe = Dsmass.frame();

	//DsmassWrong.plotOn(massframe, RooFit::Name("WrongCombination"),DrawOption("F"), FillColor(kGreen), FillStyle(2015) );
	data->plotOn(massframe, RooFit::Name("Data") );
	model->plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
	model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("bkg"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
	model->plotOn(massframe, RooFit::Name("Signal"),Components("sig"),Range("all"), LineColor(kOrange), LineStyle(kDashed), DrawOption("F"), FillColor(kOrange), FillStyle(2015));
	model->paramOn(massframe,Layout(0.75,0.99,0.75));
	massframe->getAttText()->SetTextSize(0.028);
	massframe->GetYaxis()->SetTitleOffset(1.3);
	massframe->SetXTitle("Dsmass (GeV)");

	TCanvas d;
  	d.SetTitle("");

	TPad *p1 = new TPad("p1","p1",0.0,0.27,0.82,0.99);
	p1->SetTitle("");
	p1->SetBorderMode(1); 
	p1->SetFrameBorderMode(0); 
	p1->SetBorderSize(2);

	p1->SetBottomMargin(0.10);

	p1->Draw();

	TPad *p2 = new TPad("p2","p2",0.0,0.065,0.82,0.24);
	p2->SetTitle("");
	p2->SetTopMargin(0.); 
	p2->SetBottomMargin(0.2);
   
	p2->SetBorderMode(1);
	p2->SetFrameBorderMode(0);
	p2->SetBorderSize(1); 
  
	p2->Draw();

	p1->cd();
	massframe->Draw();
	TLatex* tex11 = new TLatex(0.6,0.8," (pp) 13 TeV");
	tex11->SetNDC(kTRUE);
	tex11->SetLineWidth(2);
	tex11->SetTextSize(0.04);
	tex11->Draw();
	tex11 = new TLatex(0.6,0.85,"Preliminary Study");
	tex11->SetNDC(kTRUE);
	tex11->SetTextFont(42);
	tex11->SetTextSize(0.04);
	tex11->SetLineWidth(2);
	tex11->Draw();

	double sigma1_str = sigma1->getVal();
	double sigma1_err = sigma1->getError();
	double sigma2_str = sigma2->getVal();
	double sigma2_err = sigma2->getError();

	double chis = massframe->chiSquare();
	cout << "---------------------------------" << endl;
	cout << "chis: " << chis <<  endl;
	cout << "---------------------------------" << endl;

	/*TLatex* tex12 = new TLatex(0.15, 0.85, Form("#sigma_{exp} = %.3lf #pm %.3lf",sigma1_str,sigma1_err));
	tex12->SetNDC(kTRUE);
	tex12->SetTextFont(42);
	tex12->SetTextSize(0.04);
	tex12->Draw();*/

	TLatex* tex13 = new TLatex(0.15, 0.85, Form("#chi^{2}/DOF = %.3lf",chis));
	tex13->SetNDC(kTRUE);
	tex13->SetTextFont(42);
	tex13->SetTextSize(0.04);
	tex13->Draw();

	double totalEntries = data->sumEntries();

	//TLegend *leg = new TLegend (0.4, 0.5, 0.6, 0.7); 
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "Data", "LP");
	leg->AddEntry("Data",(Form("Entries: %2.0f", totalEntries))," ");
	leg->AddEntry("Fit","Sig + Bkg","L");
	leg->AddEntry("Signal","Signal","LP");
	leg->AddEntry("Combinatorial","Combinatorial","LP");
	leg->Draw();

	//Pull dists
	RooHist* pull_hist = massframe->pullHist("Data","Fit");
	RooPlot *pull_plot = Dsmass.frame();

	pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
	pull_plot->SetTitle("");
	pull_plot->GetXaxis()->SetTitle("");
	pull_plot->GetXaxis()->SetTitleFont(42);
	pull_plot->GetXaxis()->SetTitleSize(0.17);
	pull_plot->GetXaxis()->SetTitleOffset(1.09);
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelSize(0.15);
	pull_plot->GetXaxis()->SetLabelOffset(0.01);
	pull_plot->GetXaxis()->SetTickLength(0.13);

	pull_plot->GetYaxis()->SetTitle("Pull hist");
	pull_plot->GetYaxis()->SetTitleFont(42);  
	pull_plot->GetYaxis()->SetTitleSize(1.30);
	pull_plot->GetYaxis()->SetTitleOffset(1.09);
	pull_plot->GetYaxis()->SetLabelFont(42);
	pull_plot->GetYaxis()->SetLabelSize(0.13);
	pull_plot->GetYaxis()->SetLabelOffset(0.005);
	pull_plot->GetYaxis()->SetNdivisions(305);

	p2->cd();
	pull_plot->Draw();
	
	d.SaveAs("FIT_plot_complete.pdf");

}//End plot_complete_fit
//=====================================================================
void DoSPlot(RooWorkspace& ws){

	//we need the fit and the dataset previously saved in the woorkspace
	RooDataSet* data = (RooDataSet*) ws.data("data");
   RooAbsPdf* model = ws.pdf("model");
   
	//we need the n values previously saved in the woorkspace
   RooRealVar* D0Yield = ws.var("n_signal");
   RooRealVar* BgYield = ws.var("n_combinatorial");	
	
   //fit the model to the data
   model->fitTo(*data,Extended());

	//sPlot technique requires model parameters (other than the yields) to be fixed
   RooRealVar* mean  = ws.var("mean");
   RooRealVar* sigma1 = ws.var("sigma1");
   RooRealVar* sigma2 = ws.var("sigma2");
   RooRealVar* coef = ws.var("coef");
   RooRealVar* lambda = ws.var("lambda");	

	mean->setConstant();
   sigma1->setConstant();
   sigma2->setConstant();
   coef->setConstant();
   lambda->setConstant();
	
	//????????????????????????????????
	RooMsgService::instance().setSilentMode(true);

	//add sWeights to dataset based on model and yield variables
  //sPlot class adds a new variable that has the name of the corresponding yield + "_sw".
   SPlot* sData = new SPlot("sData","An sPlot",*data, model, RooArgList(*D0Yield,*BgYield));

	cout << endl <<  "Yield of D0 is "
       << D0Yield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_signal") << endl;

  cout << "Yield of background is "
       << BgYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_combinatorial") << endl
       << endl;

	// The first ten entries	
	for(Int_t i=0; i < 10; i++)
    {
      std::cout << "D0 Weight   " << sData->GetSWeight(i,"n_signal")
                << "   Bg Weight  " << sData->GetSWeight(i,"n_combinatorial")
                << "  Total Weight   " << sData->GetSumOfEventSWeight(i)
                << std::endl;
    }

	//the reweighted data is saved in the woorkspace
	cout << "import data with wights" << endl;
	ws.import(*data, Rename("dataWithSWeights"));


}//End DoSPlot
//=====================================================================
void MakePlots(RooWorkspace& ws){
	
	// make our canvas
	TCanvas* cdata = new TCanvas("sPlot","sPlot", 800, 600);
	cdata->Divide(2,2);

	RooAbsPdf* model = ws.pdf("model");
	RooAbsPdf* sig = ws.pdf("sig");
	RooAbsPdf* bkg = ws.pdf("bkg");

	RooRealVar* Dsmass  = ws.var("Dsmass");
	RooRealVar* D0pt = ws.var("D0pt");

	RooRealVar* BpYield = ws.var("n_signal");
	RooRealVar* BgYield = ws.var("n_combinatorial");

	double sigYield = BpYield->getVal();
	double bkgYield = BgYield->getVal();
	double n = 100;

	RooDataSet* data = (RooDataSet*) ws.data("data");

	cdata->cd(1);
	RooPlot* mframe = Dsmass->frame();
	mframe->GetXaxis()->SetTitle(TString::Format("mass of D0 [GeV]"));
	data->plotOn(mframe);
	model->plotOn(mframe,LineColor(kRed));
	model->plotOn(mframe,Components(*sig),LineStyle(kDashed),LineColor(kOrange));
	model->plotOn(mframe,Components(*bkg),LineStyle(kDashed),LineColor(kBlue));
	mframe->SetTitle("Dsmass");
	mframe->Draw();

	cdata->cd(2);
	RooPlot* ptframe = D0pt->frame();
	data->plotOn(ptframe); 
	ptframe->SetTitle("pT of D0: total sample");
	ptframe->Draw();

	//get the dataset with sWeights
	// The SPlot class adds a new variable that has the name of the corresponding yield + "_sw".
	RooDataSet* dataW = (RooDataSet*) ws.data("dataWithSWeights");
	RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
	ws.import(*dataWBp,Rename("dataWBp"));
	RooDataSet* dataWBg = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_combinatorial_sw");

	RooPlot* ptframe2Bp = D0pt->frame();
	RooPlot* ptframe2Bg = D0pt->frame();

	ptframe2Bp->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(D0pt->getMax()-D0pt->getMin())/n));
	ptframe2Bg->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(D0pt->getMax()-D0pt->getMin())/n));

	ptframe2Bp->GetXaxis()->SetTitle("pT of D0");
	ptframe2Bg->GetXaxis()->SetTitle("pT of D0");
  
  	dataWBp->plotOn(ptframe2Bp, DataError(RooAbsData::SumW2),Binning(n));
  	dataWBg->plotOn(ptframe2Bg, DataError(RooAbsData::SumW2),Binning(n));

	ptframe2Bp->SetTitle("pT distribution of D0 for signal (splot)");
	ptframe2Bg->SetTitle("pT distribution of D0 for background (splot)");
  
	cdata->cd(3);  ptframe2Bp->Draw();
	cdata->cd(4);  ptframe2Bg->Draw();

	cdata->SaveAs("FitSPlot.pdf");
  //cdata->SaveAs("SPlot.gif");


}
