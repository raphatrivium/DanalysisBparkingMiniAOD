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
#include <RooExponential.h>

using namespace RooStats;
using namespace RooFit;
using namespace std;

// see below for implementation
void set_up_workspace_variables(RooWorkspace& ws);
void AddData(RooWorkspace& ws, TString f_input);
void AddModel(RooWorkspace& ws);
void plot_complete_fit(RooWorkspace& ws);
void DoSPlot(RooWorkspace& ws);
void MakePlots(RooWorkspace& ws, int n, TString label);
void pT_analysis(RooWorkspace& ws, int n, TString ptfile, TString datafile);
void eta_analysis(RooWorkspace& ws, int n, TString ptfile, TString datafile);
void ForLifetime(RooWorkspace& ws);

//--------------------------------------
// D E F I N E  T H E  P A R T I C L E
#define particle 0 // 0= D* ; 1= D0
// D E F I N E  T H E  D A T A S E T
#define Dataset 0 // 0= Data ; 1=MC ; 2=MinBias
//--------------------------------------
void DstarD0Fit()
{	//For the out put of time taken to this program run
	clock_t tStart = clock();

	#if Dataset == 0
	TString input_file_data = "BparkingDataset_Data.root";
	#elif Dataset == 1
	TString input_file_data = "MC_DStarToD0Pi_D0KPi.root";
	#elif Dataset == 2
	TString input_file_data = "MC_MinBias.root";
	#endif

#if particle == 0
	TString variables[] = {"Dseta","Dsphi","Dspt","D0fromDSs3D",
									"D0mass","D0eta","D0phi","D0pt", "Dslifetime"};
	int n_var = sizeof(variables)/sizeof(variables[0]);
	int n_bins[]= {20, 20, 10, 10, 20, 10, 10, 10, 10, 10, 15, 10, 10, 15, 15, 
						15, 15, 15, 15, 15, 15, 20, 20, 20, 20, 20, 20, 20, 20};
#elif particle == 1
	TString variables[] = {"D0Kpipt","D0Kpieta","D0Kpiphi","D0lifetime"
									};
	int n_var = sizeof(variables)/sizeof(variables[0]);
	int n_bins[]= {20, 20, 10, 10, 20, 10, 10, 10, 10, 10, 15, 10, 10, 15, 15, 
						15, 15, 15, 15, 15, 15, 20, 20, 20, 20, 20, 20, 20, 20};
#endif
	//-------------------------------------------------------------
	// Create a new workspace to manage the project.
  	RooWorkspace* MyWS = new RooWorkspace("MyWS");
	//-------------------------------------------------------------
	// set up variables.
  	set_up_workspace_variables(*MyWS);
	//-------------------------------------------------------------
	// add some toy data to the workspace (depends of set_up_workspace_variables(*MyWS))
	AddData(*MyWS, input_file_data);
	//-------------------------------------------------------------
	// Inside this function you will find a discription our model. 
	AddModel(*MyWS);
	//-------------------------------------------------------------
	//Fit Signal + Background. (depends of AddData())
	plot_complete_fit(*MyWS);
	//-------------------------------------------------------------
	// do sPlot.This wil make a new dataset with sWeights added for every event. (depends of plot_complete_fit())
	DoSPlot(*MyWS);
	//-------------------------------------------------------------
	// Make some plots showing the discriminating variable. (depends of DoSPlot())
	for (int i=0; i < n_var; i++ ) {MakePlots(*MyWS, n_var, variables[i]);}
	//-------------------------------------------------------------
	//Inspect the workspace if you wish
	//MyWS->Print();
	//-------------------------------------------------------------
	//Eta and pT analysis (depends of plot_complete_fit())
	//pT_analysis(*MyWS,n_bins[0],"pT.root", input_file_data); 
	//eta_analysis(*MyWS,n_bins[0],"eta.root", input_file_data);
	//-------------------------------------------------------------
	//Analysis of the Lifetime (depends of DoSPlot())
	ForLifetime(*MyWS); 

  	//cleanup
  	delete MyWS;

	//T I M E
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	printf("Time taken: %.2fm\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC*60.));
	printf("Time taken: %.2fh\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC*60.*60.));		
	cout << "---------------------------------" << endl;
	cout << "E N D   P R O G R A M" << endl;
	cout << "---------------------------------" << endl;
}//end main program
//=====================================================================
void set_up_workspace_variables(RooWorkspace& ws)
{
	cout << "---------------------------------" << endl;
	cout << "set_up_workspace_variables" << endl;
	cout << "---------------------------------" << endl;

#if particle == 0
	// set range of observable
	Double_t lowRange = 1.93, highRange = 2.10;
	Double_t lowRange2 = 1.77, highRange2 = 1.95;
	Double_t lowEta = -4., highEta = 4.;
	Double_t lowPhi = -4., highPhi = 4.;
	Double_t lowdxy = -0.1, highdxy = 0.1;
	Double_t lowdz = -0.1, highdz = 0.1;
	Double_t lowlifetime = 0.*(pow(10,-12));
	Double_t highlifetime = 1.6*(pow(10,-12));

	// make a RooRealVar for the observables
	RooRealVar Dsmass("Dsmass","D* M_{inv}",lowRange,highRange, "GeV/c^{2}");
	RooRealVar Dseta("Dseta","#eta",lowEta,highEta);
	RooRealVar Dsphi("Dsphi","#phi",lowPhi,highPhi);
	RooRealVar Dspt("Dspt","Transverse Momentum",0.,100.,"GeV/c");	
	RooRealVar D0fromDSs3D("D0fromDSs3D","D0fromDSs3D",0.,5.);
	RooRealVar D0mass("D0mass","D0mass",lowRange2,highRange2);
	RooRealVar D0phi("D0phi","D0phi",lowPhi,highPhi);
	RooRealVar D0eta("D0eta","D0eta",lowEta,highEta);
	RooRealVar D0pt("D0pt","D0pt",4.,100.);
	RooRealVar Dslifetime("Dslifetime","Dslifetime",lowlifetime,highlifetime);
	#if Dataset == 0
	#elif Dataset == 1
	RooRealVar PUWeight("PUWeight","PUWeight",0.,10.,"Weight");
	ws.import(PUWeight);
	#elif Dataset == 2
	RooRealVar PUWeight("PUWeight","PUWeight",0.,10.,"Weight");
	ws.import(PUWeight);
	#endif

	RooRealVar DsmassWrong("DsmassWrong","D* M_{inv}", lowRange, highRange, "GeV/c^{2}");
	
	// import to your worksapce
	ws.import(Dsmass);
	ws.import(Dseta);
	ws.import(Dsphi);
	ws.import(Dspt);

	ws.import(D0fromDSs3D);
	ws.import(D0mass);
	ws.import(D0phi);
	ws.import(D0eta);
	ws.import(D0pt); 
	ws.import(Dslifetime);

	ws.import(DsmassWrong);

#elif particle == 1
  	// set range of observable
	Double_t lowRange = 1.77, highRange = 1.96;
	Double_t lowEta = -4., highEta = 4.;
	Double_t lowPhi = -4., highPhi = 4.;
	Double_t lowdxy = -0.1, highdxy = 0.1;
	Double_t lowdz = -0.1, highdz = 0.1;
	Double_t lowlifetime = 0.0*(pow(10,-12));
	Double_t highlifetime = 5.0*(pow(10,-12));
	
	// make a RooRealVar for the observables
	RooRealVar D0Kpimass("D0Kpimass","D0 Invariant Mass",lowRange,highRange, "GeV/c^{2}");
	RooRealVar D0Kpieta("D0Kpieta","#eta",lowEta,highEta);
	RooRealVar D0Kpiphi("D0Kpiphi","#phi",lowPhi,highPhi);
	RooRealVar D0Kpipt("D0Kpipt","Transverse Momentum",0.,100.,"GeV/c");
	RooRealVar D0lifetime("D0lifetime","D0lifetime",lowlifetime,highlifetime,"s");
	#if Dataset == 0
	#elif Dataset == 1
	RooRealVar PUWeight("PUWeight","PUWeight",0.,10.,"Weight");
	ws.import(PUWeight);
	#elif Dataset == 2
	RooRealVar PUWeight("PUWeight","PUWeight",0.,10.,"Weight");
	ws.import(PUWeight);
	#endif

	ws.import(D0Kpimass);
	ws.import(D0Kpieta);
	ws.import(D0Kpiphi);
	ws.import(D0Kpipt);
	ws.import(D0lifetime);
#endif

	cout << "---------------------------------" << endl;
	cout << "End set_up_workspace_variables" << endl;
	cout << "---------------------------------" << endl;
}
//=====================================================================
void AddData(RooWorkspace& ws, TString f_input){
	cout << "---------------------------------" << endl;
	cout << "AddData" << endl;
	cout << "---------------------------------" << endl;

	TFile* file = new TFile(f_input);

#if particle == 0
	TTree* t1_data = (TTree *)file->Get("t_analysis");
	TTree* t2_data = (TTree *)file->Get("t_DstarWrongCombination");
	RooArgList arg_list ("arg_list");

	arg_list.add(*(ws.var("Dsmass")));
	arg_list.add(*(ws.var("Dspt")));
	arg_list.add(*(ws.var("Dsphi")));
	arg_list.add(*(ws.var("Dseta")));
	arg_list.add(*(ws.var("D0fromDSs3D")));
	arg_list.add(*(ws.var("D0mass")));
	arg_list.add(*(ws.var("D0phi")));
	arg_list.add(*(ws.var("D0eta")));
	arg_list.add(*(ws.var("D0pt"))); 
	arg_list.add(*(ws.var("Dslifetime")));
	#if Dataset == 0
	#elif Dataset == 1
	RooRealVar PUWeight = *(ws.var("PUWeight"));
	#elif Dataset == 2
	RooRealVar PUWeight = *(ws.var("PUWeight"));
	#endif

	RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);

	#if Dataset == 0
	ws.import(*data, Rename("data"));
	#elif Dataset == 1
	// Add column with variable w to previously generated dataset
	RooRealVar* w = (RooRealVar*) data->addColumn(PUWeight) ;
	// Instruct dataset wdata in interpret w as event weight rather than as observable
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName()) ;
	ws.import(wdata, Rename("data"));
	#elif Dataset == 2
	RooRealVar* w = (RooRealVar*) data->addColumn(PUWeight) ;
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName()) ;
	ws.import(wdata, Rename("data"));
	#endif
	
#elif particle == 1
	TTree* t1_data = (TTree *)file->Get("t_D0analysis");
	RooArgList arg_list ("arg_list");

	arg_list.add(*(ws.var("D0Kpimass")));
	arg_list.add(*(ws.var("D0Kpieta")));
	arg_list.add(*(ws.var("D0Kpiphi")));
	arg_list.add(*(ws.var("D0Kpipt")));
	arg_list.add(*(ws.var("D0lifetime")));
	#if Dataset == 0
	#elif Dataset == 1
	RooRealVar PUWeight = *(ws.var("PUWeight"));
	#elif Dataset == 2
	RooRealVar PUWeight = *(ws.var("PUWeight"));
	#endif

	RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);
	#if Dataset == 0
	ws.import(*data, Rename("data"));
	#elif Dataset == 1
	// Add column with variable w to previously generated dataset
	RooRealVar* w = (RooRealVar*) data->addColumn(PUWeight) ;
	// Instruct dataset wdata in interpret w as event weight rather than as observable
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName()) ;
	ws.import(*data, Rename("data"));
	#elif Dataset == 2
	RooRealVar* w = (RooRealVar*) data->addColumn(PUWeight) ;
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName()) ;
	ws.import(wdata, Rename("data"));
	#endif	

#endif

	cout << "---------------------------------" << endl;
	cout << "End AddData" << endl;
	cout << "---------------------------------" << endl;
}//End addData
//=====================================================================
void AddModel(RooWorkspace& ws) {
	cout << "---------------------------------" << endl;
	cout << "AddModel" << endl;
	cout << "---------------------------------" << endl;

#if particle == 0
	RooRealVar Dsmass = *(ws.var("Dsmass"));
	RooDataSet* data = (RooDataSet*) ws.data("data");
	///------------------------------------------------
	// P D F   S I G N A L 
	// ------------------------------------------------
	RooRealVar mean( "mean", "mean", 2.010, 1.900, 2.100);
	RooRealVar sigma1( "sigma1", "sigma1", 0.0189,-2.00, 2.00);
	RooRealVar sigma2( "sigma2", "sigma2", 0.00930, -2.00, 2.00);
	RooGaussian gauss1("gauss1","gauss1",Dsmass,mean,sigma1);
	RooGaussian gauss2("gauss2","gauss2",Dsmass,mean,sigma2);
	RooRealVar coef("coef","fraction of Gaussian1", 0.452, 0.0, 1.0);
	RooAddPdf  sig("sig","gauss1 + gauss2",RooArgList(gauss1,gauss2),RooArgList(coef));
	//Make the RooRealVar mean constant. It will not vary in the input range
	//mean.setConstant();	

	#if Dataset == 0
	sigma1.setConstant();
	sigma2.setConstant();
	#elif Dataset == 1
	#elif Dataset == 2
	#endif
	
	

	//RooRealVar tail( "tail", "tail", 0.00912, -10.00, 10.00);
	//RooNovosibirsk novosibirsk("novosibirsk", "novosibirsk", Dsmass, mean, sigma1, tail);
	//RooAddPdf  sig("sig","novosibirsk1",RooArgList(novosibirsk));	
	//------------------------------------------------
	// P D F   B A C K G R O U N D
	// ------------------------------------------------
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

#elif particle == 1
	RooRealVar D0Kpimass = *(ws.var("D0Kpimass"));
	RooDataSet* data = (RooDataSet*) ws.data("data");
	///------------------------------------------------
	// P D F   S I G N A L 
	// ------------------------------------------------
	RooRealVar mean( "mean", "mean", 1.864, 1.840, 1.900);
	RooRealVar sigma1( "sigma1", "sigma1", 0.01740,-2.00, 2.00);
	RooRealVar sigma2( "sigma2", "sigma2", 0.00912, -2.00, 2.00);
	RooGaussian gauss1("gauss1","gauss1",D0Kpimass,mean,sigma1);
	RooGaussian gauss2("gauss2","gauss2",D0Kpimass,mean,sigma2);
	RooRealVar coef("coef","fraction of Gaussian1", 0.452, 0.0, 1.0);
	RooAddPdf  sig("sig","gauss1 + gauss2",RooArgList(gauss1,gauss2),RooArgList(coef));
	// Make the RooRealVar mean constant. It will not vary in the input range
	//mean.setConstant();
   //sigma1.setConstant();
   //sigma2.setConstant();		
	//------------------------------------------------
	// P D F   B A C K G R O U N D
	// ------------------------------------------------
	RooRealVar lambda("lambda","lambda",-2.26,-10.,10.0);
  	RooExponential bkg("bkg", "bkg", D0Kpimass, lambda);
	/*RooRealVar poly_c1("poly_c1", "coefficient of x^1 term", 0, -10, 10);
	RooRealVar poly_c2("poly_c2", "coefficient of x^2 term", 0, -10, 10);
	RooRealVar poly_c3("poly_c3", "coefficient of x^3 term", 0, -10, 10);
	RooPolynomial bkg("bkg", "bkg", D0Kpimass, RooArgList(poly_c1, poly_c2, poly_c3) );*/
	//------------------------------------------------
	// A d d  s i g n a l   a n d   b a c k g r o u n d
	// ------------------------------------------------
	double n_signal_initial = data->sumEntries(TString::Format("abs(D0Kpimass-%g)<0.05",mean.getVal())) - data->sumEntries(TString::Format("abs(D0Kpimass-%g)<0.10&&abs(D0Kpimass-%g)>0.05",mean.getVal(),mean.getVal()));
	double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
	RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
	RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
	RooAddPdf  model("model","sig + bkg", RooArgList(sig,bkg), RooArgList(n_signal,n_combinatorial));
#endif
	ws.import(model);

	cout << "---------------------------------" << endl;
	cout << "End AddModel" << endl;
	cout << "---------------------------------" << endl;
}//End AddModel
//=====================================================================
void plot_complete_fit(RooWorkspace& ws){
	cout << "---------------------------------" << endl;
	cout << "plot_complete_fit" << endl;
	cout << "---------------------------------" << endl;

	RooAbsPdf* model = ws.pdf("model");
	RooDataSet* data = (RooDataSet*) ws.data("data");
	RooRealVar* sigma1 = ws.var("sigma1");
	RooRealVar* sigma2 = ws.var("sigma2");
#if particle == 0
	RooRealVar Dsmass = *(ws.var("Dsmass"));
	RooPlot* massframe = Dsmass.frame();
	massframe->SetTitle("D* Invariant Mass");
#elif particle == 1
	RooRealVar D0Kpimass = *(ws.var("D0Kpimass"));
	RooPlot* massframe = D0Kpimass.frame();
	massframe->SetTitle("D0 Invariant Mass");
#endif
	//-------------------------------
	//F I T
	//-------------------------------
	model->fitTo(*data,Range("all"));
	//-------------------------------
	//P L O T
	//-------------------------------
	//DsmassWrong.plotOn(massframe, RooFit::Name("WrongCombination"),DrawOption("F"), FillColor(kGreen), FillStyle(2015) );
	data->plotOn(massframe, RooFit::Name("Data") );
	model->plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
	model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("bkg"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
	model->plotOn(massframe, RooFit::Name("Signal"),Components("sig"),Range("all"), LineColor(kOrange), LineStyle(kDashed), DrawOption("F"), FillColor(kOrange), FillStyle(2015));
	model->paramOn(massframe,Layout(0.75,0.99,0.75));
	massframe->getAttText()->SetTextSize(0.028);
	massframe->GetYaxis()->SetTitleOffset(1.3);
	massframe->SetXTitle("Invariant Mass [GeV/c^{2}]");


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
	#if Dataset == 0
	TLatex* tex11 = new TLatex(0.6,0.8,"Data (pp) 13 TeV");
	#elif Dataset == 1
	TLatex* tex11 = new TLatex(0.6,0.8,"MC (pp) 13 TeV");
	#elif Dataset == 2
	TLatex* tex11 = new TLatex(0.6,0.8,"MinBias (pp) 13 TeV");
	#endif
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

	/*TLatex* tex12 = new TLatex(0.15, 0.85, Form("#sigma_{exp} = %.3lf #pm %.3lf",sigma1_str,sigma1_err));
	tex12->SetNDC(kTRUE); tex12->SetTextFont(42); tex12->SetTextSize(0.04);
	tex12->Draw();*/

	TLatex* tex13 = new TLatex(0.15, 0.85, Form("#chi^{2}/DOF = %.3lf",chis));
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04);
	tex13->Draw();

	double totalEntries = data->sumEntries();

#if particle == 0
	#if Dataset == 0
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "Data", "LP");
	#elif Dataset == 1
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "MC", "LP");
	#elif Dataset == 2
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "MinBias", "LP");
	#endif
#elif particle == 1
	#if Dataset == 0
	TLegend *leg = new TLegend (0.15,0.25,0.35,0.45);
	leg->AddEntry(massframe->findObject("Data"), "Data", "LP");
	#elif Dataset == 1
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "MC", "LP");
	#elif Dataset == 2
	TLegend *leg = new TLegend (0.15,0.25,0.35,0.45);
	leg->AddEntry(massframe->findObject("Data"), "MinBias", "LP");
	#endif		
#endif

	//TLegend *leg = new TLegend (0.4, 0.5, 0.6, 0.7); 
	#if Dataset == 0
	
	#elif Dataset == 1
	
	#elif Dataset == 2
	
	#endif
	leg->AddEntry("Data",(Form("Entries: %2.0f", totalEntries))," ");
	leg->AddEntry("Fit","Sig + Bkg","L");
	leg->AddEntry("Signal","Signal","LP");
	leg->AddEntry("Combinatorial","Combinatorial","LP");
	leg->Draw();

	//Pull dists
	RooHist* pull_hist = massframe->pullHist("Data","Fit");
#if particle == 0
	RooPlot *pull_plot = Dsmass.frame();
#elif particle == 1
	RooPlot *pull_plot = D0Kpimass.frame();
#endif


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

#if particle == 0
	#if Dataset == 0
		d.SaveAs("DstarData/FitDstarDATA.pdf");
	#elif Dataset == 1
		d.SaveAs("DstarMC/FitDstarMC.pdf");
	#elif Dataset == 2
		d.SaveAs("DstarMinBias/FitDstarMinBias.pdf");
	#endif
#elif particle == 1
	#if Dataset == 0
		d.SaveAs("D0DATA/FitD0DATA.pdf");
	#elif Dataset == 1
		d.SaveAs("D0MC/FitD0MC.pdf");
	#elif Dataset == 2
		d.SaveAs("D0MinBias/FitD0MinBias.pdf");
	#endif		
#endif

	cout << "---------------------------------" << endl;
	cout << "End plot_complete_fit" << endl;
	cout << "---------------------------------" << endl;
}//End plot_complete_fit
//=====================================================================
void DoSPlot(RooWorkspace& ws){
	cout << "---------------------------------" << endl;
	cout << "DoSPlot" << endl;
	cout << "---------------------------------" << endl;

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

	cout << endl <<  "Yield of Ds is "
       << D0Yield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_signal") << endl;

  cout << "Yield of background is "
       << BgYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_combinatorial") << endl
       << endl;

	// The first ten entries	
	for(Int_t i=0; i < 10; i++)
	{
      std::cout << "Ds Weight   " << sData->GetSWeight(i,"n_signal")
                << "   Bg Weight  " << sData->GetSWeight(i,"n_combinatorial")
                << "  Total Weight   " << sData->GetSumOfEventSWeight(i)
                << std::endl;
	}

	//the reweighted data is saved in the woorkspace
	ws.import(*data, Rename("dataWithSWeights"));

	cout << "---------------------------------" << endl;
	cout << "End DoSPlot" << endl;
	cout << "---------------------------------" << endl;
}//End DoSPlot
//=====================================================================
void MakePlots(RooWorkspace& ws, int n, TString label){

	// make our canvas
	TCanvas* cdata = new TCanvas("sPlot","sPlot", 800, 600);
	cdata->Divide(2,2);

	RooAbsPdf* model = ws.pdf("model");
	RooAbsPdf* sig = ws.pdf("sig");
	RooAbsPdf* bkg = ws.pdf("bkg");
	RooRealVar* variable = ws.var(label);

	#if particle == 0
	RooRealVar* Dsmass  = ws.var("Dsmass");
	#elif particle == 1
	RooRealVar* D0Kpimass  = ws.var("D0Kpimass");
	#endif	

	RooRealVar* BpYield = ws.var("n_signal");
	RooRealVar* BgYield = ws.var("n_combinatorial");

	double sigYield = BpYield->getVal();
	double bkgYield = BgYield->getVal();
	double nbin = 100;

	RooDataSet* data = (RooDataSet*) ws.data("data");

	cdata->cd(1);
#if particle == 0
	RooPlot* mframe = Dsmass->frame();
	mframe->SetTitle("D* Mass Distribuition");
#elif particle == 1
	RooPlot* mframe = D0Kpimass->frame();
	mframe->SetTitle("D0 Mass Distribuition");
#endif	

	mframe->GetXaxis()->SetTitle(TString::Format("Invariant Mass [GeV/c^{2}]"));
	data->plotOn(mframe);
	model->plotOn(mframe,LineColor(kRed));
	model->plotOn(mframe,Components(*sig),LineStyle(kDashed),LineColor(kOrange));
	model->plotOn(mframe,Components(*bkg),LineStyle(kDashed),LineColor(kBlue));
	mframe->Draw();

	cdata->cd(2);
	RooPlot* ptframe = variable->frame();
	data->plotOn(ptframe);
	#if Dataset == 0
	ptframe->SetTitle(label + ": Total Sample (Data)");
	#elif Dataset == 1
	ptframe->SetTitle(label + ": Total Sample (MC)");
	#elif Dataset == 2
	ptframe->SetTitle(label + ": Total Sample (MinBias)");
	#endif 

	ptframe->Draw();

	//get the dataset with sWeights
	// The SPlot class adds a new variable that has the name of the corresponding yield + "_sw".
	RooDataSet* dataW = (RooDataSet*) ws.data("dataWithSWeights");
	RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
	RooDataSet* dataWBg = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_combinatorial_sw");

	RooPlot* ptframe2Bp = variable->frame();
	RooPlot* ptframe2Bg = variable->frame();

	ptframe2Bp->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/nbin));
	ptframe2Bg->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/nbin));

	ptframe2Bp->GetXaxis()->SetTitle(label);
	ptframe2Bg->GetXaxis()->SetTitle(label);
  
  	dataWBp->plotOn(ptframe2Bp, DataError(RooAbsData::SumW2),Binning(nbin));
  	dataWBg->plotOn(ptframe2Bg, DataError(RooAbsData::SumW2),Binning(nbin));

	#if Dataset == 0
	ptframe2Bp->SetTitle(label + " Signal-Data (splot)");
	ptframe2Bg->SetTitle(label + " Background-Data (splot)");
	#elif Dataset == 1
	ptframe2Bp->SetTitle(label + " Signal-MC (splot)");
	ptframe2Bg->SetTitle(label + " Background-MC (splot)");
	#elif Dataset == 2
	ptframe2Bp->SetTitle(label + " Signal-MinBias (splot)");
	ptframe2Bg->SetTitle(label + " Background-MinBias (splot)");
	#endif 
	 
	cdata->cd(3);  ptframe2Bp->Draw();
	cdata->cd(4);  ptframe2Bg->Draw();

#if particle == 0
	#if Dataset == 0
	cdata->SaveAs("DstarData/SPlot"+label+"DATA.pdf");
	#elif Dataset == 1
	cdata->SaveAs("DstarMC/SPlot"+label+"MC.pdf");
	#elif Dataset == 2
	cdata->SaveAs("DstarMinBias/SPlot"+label+"MinBias.pdf");
	#endif
#elif particle == 1
	#if Dataset == 0
	cdata->SaveAs("D0DATA/SPlot"+label+"DATA.pdf");
	#elif Dataset == 1
	cdata->SaveAs("D0MC/SPlot"+label+"MC.pdf");
	#elif Dataset == 2
	cdata->SaveAs("D0MinBias/SPlot"+label+"MinBias.pdf");
	#endif		
#endif

	cout << "---------------------------------" << endl;
	cout << "End MakePlots" << endl;
	cout << "---------------------------------" << endl;
}
//=====================================================================
void build_pdf(RooWorkspace& w, std::string choice) {

#if particle == 0
	RooRealVar Dsmass = *(w.var("Dsmass"));
#elif particle == 1
	RooRealVar D0Kpimass = *(w.var("D0Kpimass"));	
#endif

  RooDataSet* data = (RooDataSet*) w.data("data");
  //  RooDataSet* reduceddata_central;
  double left =  5.3 ;
  double right = 5.45 ;
  double mass_peak = 2.010 ;

	//  reduceddata_central = (RooDataSet*) data->reduce(Form("Dsmass>%lf",left));
	//  reduceddata_central = (RooDataSet*) reduceddata_central->reduce(Form("Dsmass<%lf",right));

	//----------------------------------
	//S I G N A L
	//----------------------------------
	RooRealVar mean("mean","mean",mass_peak,mass_peak-0.1,mass_peak+0.1);
	RooRealVar sigma1("sigma1","sigma1",0.0252,0.020,0.030);
	RooRealVar sigma2("sigma2","sigma2",0.01052,0.010,0.020);
#if particle == 0
	RooGaussian signal1("signal1","signal_gauss1",Dsmass,mean,sigma1);
	RooGaussian signal2("signal2","signal_gauss2",Dsmass,mean,sigma2);
#elif particle == 1
	RooGaussian signal1("signal1","signal_gauss1",D0Kpimass,mean,sigma1);
	RooGaussian signal2("signal2","signal_gauss2",D0Kpimass,mean,sigma2);	
#endif
  RooRealVar cofs("cofs", "cofs", 0.317, 0., 1.);
  RooAddPdf signal("signal", "signal", RooArgList(signal1,signal2),cofs);
  //sigma1.setConstant();
  //sigma2.setConstant();
  //cofs.setConstant();

	//---------------------------
	//systematics->signal
	//--------------------------
	RooRealVar sigma3("sigma3","sigma3",0.012,0.010,0.030);
#if particle == 0
	RooGaussian signal3("signal3","signal3",Dsmass, mean, sigma3);
#elif particle == 1
	RooGaussian signal3("signal3","signal3",D0Kpimass, mean, sigma3);
#endif
	//----------------------------
	//BACKGROUND//
	//---------------------------------------------
	//error function
	RooRealVar m_nonprompt_scale("m_nonprompt_scale", "m_nonprompt_scale", 4.74168e-02, 0, 1);
	RooRealVar m_nonprompt_shift("m_nonprompt_shift", "m_nonprompt_shift", 5.14425, 4.5, 6.);
  
	m_nonprompt_shift.setConstant(kTRUE);
	m_nonprompt_scale.setConstant(kTRUE);

#if particle == 0
	RooGenericPdf erf("erf","erf", "TMath::Erfc((Dsmass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(Dsmass, m_nonprompt_scale, m_nonprompt_shift));
#elif particle == 1
	RooGenericPdf erf("erf","erf", "TMath::Erfc((D0Kpimass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(D0Kpimass, m_nonprompt_scale, m_nonprompt_shift));
#endif

	//exponential
	RooRealVar lambda("lambda","lambda",-2.,-5.,0.0);
#if particle == 0
	RooExponential fit_side("fit_side", "fit_side_exp", Dsmass, lambda);
#elif particle == 1
	RooExponential fit_side("fit_side", "fit_side_exp", D0Kpimass, lambda);
#endif

  //jpsi_pi component
#if particle == 0
  RooRealVar m_jpsipi_mean1("m_jpsipi_mean1","m_jpsipi_mean1",5.34693e+00,Dsmass.getAsymErrorLo(),Dsmass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean2("m_jpsipi_mean2","m_jpsipi_mean2",5.46876e+00,Dsmass.getAsymErrorLo(),Dsmass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean3("m_jpsipi_mean3","m_jpsipi_mean3",5.48073e+00,Dsmass.getAsymErrorLo(),Dsmass.getAsymErrorHi());
#elif particle == 1
  RooRealVar m_jpsipi_mean1("m_jpsipi_mean1","m_jpsipi_mean1",5.34693e+00,D0Kpimass.getAsymErrorLo(),D0Kpimass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean2("m_jpsipi_mean2","m_jpsipi_mean2",5.46876e+00,D0Kpimass.getAsymErrorLo(),D0Kpimass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean3("m_jpsipi_mean3","m_jpsipi_mean3",5.48073e+00,D0Kpimass.getAsymErrorLo(),D0Kpimass.getAsymErrorHi());
#endif

  RooRealVar m_jpsipi_sigma1l("m_jpsipi_sigma1l","m_jpsipi_sigma1l",2.90762e-02,0.010,0.150);
  RooRealVar m_jpsipi_sigma1r("m_jpsipi_sigma1r","m_jpsipi_sigma1r",6.52519e-02,0.010,0.150);
  RooRealVar m_jpsipi_sigma2("m_jpsipi_sigma2","m_jpsipi_sigma2",9.94712e-02,0.020,0.500);

  RooRealVar m_jpsipi_sigma3("m_jpsipi_sigma3","m_jpsipi_sigma3",3.30152e-01,0.020,0.500);
  RooRealVar m_jpsipi_fraction2("m_jpsipi_fraction2","m_jpsipi_fraction2",2.34646e-01,0.0,1.0);
  RooRealVar m_jpsipi_fraction3("m_jpsipi_fraction3","m_jpsipi_fraction3",1.14338e-01,0.0,1.0);

	m_jpsipi_mean1.setConstant(kTRUE);
	m_jpsipi_mean2.setConstant(kTRUE);
	m_jpsipi_mean3.setConstant(kTRUE);
	m_jpsipi_sigma1l.setConstant(kTRUE);
	m_jpsipi_sigma1r.setConstant(kTRUE);
	m_jpsipi_sigma2.setConstant(kTRUE);
	m_jpsipi_sigma3.setConstant(kTRUE);
	m_jpsipi_fraction2.setConstant(kTRUE);
	m_jpsipi_fraction3.setConstant(kTRUE);

#if particle == 0
	RooBifurGauss m_jpsipi_gaussian1("m_jpsipi_gaussian1","m_jpsipi_gaussian1",Dsmass,m_jpsipi_mean1,m_jpsipi_sigma1l,m_jpsipi_sigma1r);
	RooGaussian m_jpsipi_gaussian2("m_jpsipi_gaussian2","m_jpsipi_gaussian2",Dsmass,m_jpsipi_mean2,m_jpsipi_sigma2);
	RooGaussian m_jpsipi_gaussian3("m_jpsipi_gaussian3","m_jpsipi_gaussian3",Dsmass,m_jpsipi_mean3,m_jpsipi_sigma3);
#elif particle == 1
  	RooBifurGauss m_jpsipi_gaussian1("m_jpsipi_gaussian1","m_jpsipi_gaussian1",D0Kpimass,m_jpsipi_mean1,m_jpsipi_sigma1l,m_jpsipi_sigma1r);
	RooGaussian m_jpsipi_gaussian2("m_jpsipi_gaussian2","m_jpsipi_gaussian2",D0Kpimass,m_jpsipi_mean2,m_jpsipi_sigma2);
	RooGaussian m_jpsipi_gaussian3("m_jpsipi_gaussian3","m_jpsipi_gaussian3",D0Kpimass,m_jpsipi_mean3,m_jpsipi_sigma3);
#endif

	RooAddPdf jpsipi("jpsipi","jpsipi",RooArgList(m_jpsipi_gaussian3,m_jpsipi_gaussian2,
	m_jpsipi_gaussian1),RooArgList(m_jpsipi_fraction3,m_jpsipi_fraction2));

#if particle == 0
	Dsmass.setRange("all", Dsmass.getMin(),Dsmass.getMax());
	Dsmass.setRange("right",right,Dsmass.getMax());
	Dsmass.setRange("left",Dsmass.getMin(),left);
	Dsmass.setRange("peak",left,right);
	Dsmass.setRange("peakright",left,Dsmass.getMax());
//n values
  double n_signal_initial = data->sumEntries(TString::Format("abs(Dsmass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(Dsmass-%g)<0.10&&abs(Dsmass-%g)>0.05",mass_peak,mass_peak));
#elif particle == 1
  	D0Kpimass.setRange("all", D0Kpimass.getMin(),D0Kpimass.getMax());
	D0Kpimass.setRange("right",right,D0Kpimass.getMax());
	D0Kpimass.setRange("left",D0Kpimass.getMin(),left);
	D0Kpimass.setRange("peak",left,right);
	D0Kpimass.setRange("peakright",left,D0Kpimass.getMax());
//n values
  double n_signal_initial = data->sumEntries(TString::Format("abs(D0Kpimass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(D0Kpimass-%g)<0.10&&abs(D0Kpimass-%g)>0.05",mass_peak,mass_peak));
#endif


  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
			cout << "build_pdf 14" << endl;
  RooRealVar f_erf("f_erf","f_erf",2.50259e-01,0,1);
  RooProduct n_erf("n_erf","n_erf",RooArgList(n_signal,f_erf));
  
  RooRealVar f_jpsipi("f_jpsipi","f_jpsipi",4.1E-5/1.026E-3,0.,0.1); 
  f_jpsipi.setConstant(kTRUE);
  RooProduct n_jpsipi("n_jpsipi","n_jpsipi",RooArgList(n_signal,f_jpsipi));
			cout << "build_pdf 15" << endl;
  //systematics->bkg
  RooRealVar slope("slope","slope",0,-10,10);
#if particle == 0
	  RooPolynomial poly_bkg("poly_bkg", "poly_bkg", Dsmass, slope);
#elif particle == 1
  	  RooPolynomial poly_bkg("poly_bkg", "poly_bkg", D0Kpimass, slope);
#endif

			cout << "build_pdf 16" << endl;
    if (choice == "nominal"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    }else if (choice=="bkg_poly"){
      RooAddPdf model("model", "model", RooArgList(signal,poly_bkg,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    } else if (choice=="bkg_range"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side, jpsipi),RooArgList(n_signal,n_combinatorial, n_jpsipi));
      w.import(model);
    }else if (choice == "signal1gauss"){
      RooAddPdf model("model", "model", RooArgList(signal3,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
      w.import(model);
    }

 			cout << "build_pdf 17" << endl;
}
//=====================================================================
double get_yield_syst(RooDataSet* data_bin, TString syst_src) {
  //returns the yield's value per bin
  
	//cout << "aaa 0\n";
	//data_bin->Print();
	//cout << "aaa 1\n";

	RooWorkspace* ws = new RooWorkspace("ws");
	set_up_workspace_variables(*ws);
	ws->import(*data_bin, Rename("data"));

	//  TString mdl = (syst_src.Contains("range")) ? "nominal" : syst_src;
	TString rng = (syst_src.Contains("range")) ? "peakright" : "all";
	//cout << "bla 0\end ";
	build_pdf(*ws,syst_src.Data());
	
	RooAbsPdf* model = ws->pdf("model");
	RooFitResult* fitres = model->fitTo(*data_bin,Range(rng),Save());
	cout << "Debug get_yield_syst 4" << endl;
	//build_pdf(*ws,"bkg_range");
	//model = ws->pdf("model");
	//RooFitResult* fitres_bgrange = model->fitTo(*data_bin,Range("peakright"),Save());

	RooRealVar* n1_var = (RooRealVar*) fitres ->floatParsFinal().find("n_signal");

	double n1  = n1_var->getVal();

	return n1; 
		
}
//=====================================================================
void pT_analysis(RooWorkspace& ws, int n, TString ptfile, TString datafile){

	TString dir_name = ".";
	TFile* f_wei = new TFile(dir_name + ptfile, "recreate");
	//dosp = true -> the splot technique is applied
	bool dosp = true; // true  -  false

	RooAbsPdf* model = ws.pdf("model");
  	RooDataSet* data = (RooDataSet*) ws.data("data");
#if particle == 0
	RooRealVar* Dspt = ws.var("Dspt");
  	RooRealVar Dsmass = *(ws.var("Dsmass"));
#elif particle == 1
	RooRealVar* D0Kpipt = ws.var("D0Kpipt");
  	RooRealVar D0Kpimass = *(ws.var("D0Kpimass"));
#endif
  	
	double mean_w = 0.;
	double mean_s = 0.;

	const int n_pt_bins = 9;
  	double pt_bins [n_pt_bins + 1] = {4,5,6,7,8,12,16,24,40,100};

	double pt_mean[n_pt_bins];
	double pt_low[n_pt_bins];
	double pt_high[n_pt_bins];

	double yield[n_pt_bins];
	double yield_err_low[n_pt_bins];
	double yield_err_high[n_pt_bins];
	double yield_err_syst[n_pt_bins];
	double m_yield_err_syst[n_pt_bins];

	for(int k = 1; k<n_pt_bins; k++){ yield_err_syst[k] = 0; m_yield_err_syst[k] = 0; }

	RooDataSet* data_pt, data_w, data_wp;
	RooFitResult* fit_pt;
	RooRealVar* n_sig_pt;
	RooRealVar* n_comb_pt;

	const int n_pdf_syst=4;
	//number of pdfs
	TString syst_src[n_pdf_syst]={"nominal","bkg_poly","signal1gauss","bkg_range"};
	//array of pdfs to evaluate the systematic uncertainty of the N signal
	double yield_syst[n_pt_bins][n_pdf_syst]; 
	//value of the systematic uncertainty

	for(int i=0;i<n_pt_bins;i++)
	{
		int int1 = (int)pt_bins[i];
		int int2 = (int)pt_bins[i+1];

		//select data subset corresponding to pT bin
#if particle == 0
		data_pt = (RooDataSet*) data->reduce(Form("Dspt>%lf",pt_bins[i]));
		data_pt = (RooDataSet*) data_pt->reduce(Form("Dspt<%lf",pt_bins[i+1]));
#elif particle == 1
		data_pt = (RooDataSet*) data->reduce(Form("D0Kpipt>%lf",pt_bins[i]));
		data_pt = (RooDataSet*) data_pt->reduce(Form("D0Kpipt<%lf",pt_bins[i+1]));
#endif
			
		ws.import(*data_pt, Rename(Form("data_pt_%d",i)));
   
		//perform fit and save result
		fit_pt = model->fitTo(*data_pt, Minos(true), Save());

		//plots the fit result
#if particle == 0
		RooPlot* massframe = Dsmass.frame(Title(""));
#elif particle == 1
		RooPlot* massframe = D0Kpimass.frame(Title(""));
#endif
		data_pt->plotOn(massframe);
		model->paramOn(massframe,Layout(0.60,0.99,0.89));
		model->plotOn(massframe, Range("all"));
	
		TCanvas b;
		massframe->Draw();
		TLatex* tex11 = new TLatex(0.15,0.75, Form("%i<pt<%i GeV",int1,int2));
		tex11->SetNDC(kTRUE); tex11->SetLineWidth(2);
		tex11->SetTextSize(0.04); tex11->Draw();
		tex11 = new TLatex(0.15,0.8," (pp) 13 TeV");
		tex11->SetNDC(kTRUE); tex11->SetTextFont(42);
		tex11->SetLineWidth(2);	tex11->SetTextSize(0.04);
		tex11->Draw();
		tex11 = new TLatex(0.15,0.85,"Preliminary Study");
		tex11->SetNDC(kTRUE); tex11->SetTextFont(42);
		tex11->SetTextSize(0.04); tex11->SetLineWidth(2);
		tex11->Draw();

#if particle == 0
	#if Dataset == 0
		b.SaveAs(Form("DstarData/PT%i_%iDstarDATA.pdf", int1, int2));
		b.SaveAs(Form("DstarData/PT%i_%iDstarDATA.png", int1, int2));
	#elif Dataset == 1
		b.SaveAs(Form("DstarMC/PT%i_%iDstarMC.pdf", int1, int2));
		b.SaveAs(Form("DstarMC/PT%i_%iDstarMC.png", int1, int2));
	#elif Dataset == 2
		b.SaveAs(Form("DstarMinBias/PT%i_%iDstarMinBias.pdf", int1, int2));
		b.SaveAs(Form("DstarMinBias/PT%i_%iDstarMinBias.png", int1, int2));
	#endif
#elif particle == 1
	#if Dataset == 0
		b.SaveAs(Form("D0DATA/PT%i_%iD0DATA.pdf", int1, int2));
		b.SaveAs(Form("D0DATA/PT%i_%iD0DATA.png", int1, int2));
	#elif Dataset == 1
		b.SaveAs(Form("D0MC/PT%i_%iD0MC.pdf", int1, int2));
		b.SaveAs(Form("D0MC/PT%i_%iD0MC.png", int1, int2));
	#elif Dataset == 2
		b.SaveAs(Form("D0MinBias/PT%i_%iD0MinBias.pdf", int1, int2));
		b.SaveAs(Form("D0MinBias/PT%i_%iD0MinBias.png", int1, int2));
	#endif		
#endif

		//get yield and its errors
		//signal yield
		n_sig_pt = (RooRealVar*) fit_pt->floatParsFinal().find("n_signal");
		//combinatorial background yield
		n_comb_pt = (RooRealVar*) fit_pt->floatParsFinal().find("n_combinatorial");

		yield[i] = n_sig_pt->getVal();
		yield_err_low [i] = n_sig_pt->getError(); 
		yield_err_high[i] = n_sig_pt->getError(); 

		cout << "test asym error:" << n_sig_pt->getErrorLo() << " " <<  n_sig_pt->getAsymErrorLo() << " symmetric: " <<  n_sig_pt->getError() <<  endl;
		
		for(int k = 1; k<n_pdf_syst; k++)
		{
			double val = 0.; 
			double val_nominal = 0.;
			val = get_yield_syst(data_pt, syst_src[k]);
			val_nominal = get_yield_syst(data_pt, syst_src[0]);
			cout << "syst nominal: "<< syst_src[0] <<endl;
			yield_syst[i][k] = (val - val_nominal);
			cout << "bin:" << i << " range bin min" << pt_bins[i] << " src:" << k << " syst:" << syst_src[k] << " yield syst:" << yield_syst[i][k]<< " rel:"<<  ((val - val_nominal)/val_nominal)*100 << "\%"<<endl;
		}  

    	for(int k = 1; k<n_pdf_syst; k++){ yield_err_syst[i] += pow(yield_syst[i][k],2); }

		m_yield_err_syst[i] = sqrt(yield_err_syst[i]);

		if (dosp)
		{
		   RooRealVar* mean  = ws.var("mean");
			RooRealVar* sigma1 = ws.var("sigma1");
			RooRealVar* sigma2 = ws.var("sigma2");
			RooRealVar* coef = ws.var("coef");
			RooRealVar* lambda = ws.var("lambda");	

			//sPlot technique requires model parameters (other than the yields) to be fixed
			mean->setConstant();
			sigma1->setConstant();
			sigma2->setConstant();
			coef->setConstant();
			lambda->setConstant();
		   
		   SPlot("sData","An sPlot",*data_pt, model, RooArgList(*n_sig_pt,*n_comb_pt));
		   
		   ws.import(*data_pt, Rename(Form("data_pt_WithSWeights_%d",i)));
		   
		   RooDataSet* data_w = (RooDataSet*) ws.data(Form("data_pt_WithSWeights_%d",i));
		   RooDataSet* data_wb = new RooDataSet(data_w->GetName(),data_w->GetTitle(),data_w,*data_w->get(),0,"n_signal_sw");
		   
		   //weighted average pT
#if particle == 0
			mean_w = data_wb->mean(*Dspt);
		   mean_s = data_pt->mean(*Dspt);
		   pt_mean[i] = data_wb->mean(*Dspt);
#elif particle == 1
			mean_w = data_wb->mean(*D0Kpipt);
		   mean_s = data_pt->mean(*D0Kpipt);
		   pt_mean[i] = data_wb->mean(*D0Kpipt);
#endif	   
		   cout<<"mean_weight:"<< mean_w <<endl;
		   cout<<"mean:"<< mean_s << endl;
		   
		   pt_low[i]= pt_mean[i]-pt_bins[i];
		   pt_high[i]= pt_bins[i+1]-pt_mean[i];
		} 
		else
		{
	   	pt_low [i]= 0.5*(pt_bins[i+1]-pt_bins[i]);
	   	pt_high[i]= 0.5*(pt_bins[i+1]-pt_bins[i]);  
		}
    
		//normalize yield to bin width
		double bin_width = pt_bins[i+1]-pt_bins[i];
		yield[i] = yield[i]/bin_width;
		yield_err_low[i] = yield_err_low[i]/bin_width;
		yield_err_high[i] = yield_err_high[i]/bin_width;
		m_yield_err_syst[i] = m_yield_err_syst[i]/bin_width;

		cout<<"pt: "<< pt_bins[i]<<"-" << pt_bins[i+1] << " mean_weight:"<< mean_w	<< "  Nsig:" <<yield[i]<< "-" << yield_err_low[i] <<"+" << yield_err_high[i]<<endl;
		cout <<"pt: "<< pt_bins[i]<<"-" << pt_bins[i+1] <<" systematic uncertainty" << m_yield_err_syst[i]<<endl;

	}

	//plot yield vs average pT
	TCanvas c;
	TMultiGraph* mg = new TMultiGraph();

	TGraphAsymmErrors* gr = new TGraphAsymmErrors(n_pt_bins,pt_mean,yield,pt_low,pt_high,yield_err_low,yield_err_high);
	gr->SetTitle("");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(1);
	gr->SetLineColor(1);
	gr->GetXaxis()->SetTitle("p_{T}(D*) [GeV/c]");
	gr->GetYaxis()->SetTitle("raw yield [GeV^{-1}]");
	//gr->Draw("AP");
	gr->Write();
	//f_wei->Close();
	//delete f_wei;

	double pt_zero[n_pt_bins];
	for (int i=0;i<n_pt_bins;i++) pt_zero[i]= 0.;

	TGraphAsymmErrors* grs = new TGraphAsymmErrors(n_pt_bins,pt_mean,yield,pt_zero,pt_zero, m_yield_err_syst, m_yield_err_syst);
	grs->SetTitle("");
	grs->SetMarkerColor(4);
	grs->SetMarkerStyle(1);
	grs->SetLineColor(2);
	grs->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
	grs->GetYaxis()->SetTitle("raw yield [GeV^{-1}]");
	grs->Write();
	f_wei->Close();
	//grs->Draw("same");

   mg->Add(gr);
   mg->Add(grs);
   mg->Draw("AP");
   mg->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
   mg->GetYaxis()->SetTitle("raw yield [GeV^{-1}]");
   
  	if (dosp) c.SaveAs("DstarData/raw_yield_pt_Ds.pdf");
	if (dosp) c.SaveAs("DstarData/raw_yield_pt_Ds.png");

	TCanvas l;
	//log scale
	l.SetLogx();
	l.SetLogy();
	TGraphAsymmErrors* grlog = new TGraphAsymmErrors(n_pt_bins,pt_mean,yield,pt_low,pt_high,yield_err_low,yield_err_high);
	grlog->SetTitle("");
	grlog->SetMarkerColor(4);
	grlog->SetMarkerStyle(21);
	grlog->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
	grlog->GetYaxis()->SetTitle("raw yield [GeV^{-1}]");
	grlog->Draw("AP");

	if (dosp) l.SaveAs("DstarData/raw_yield_pt_logscale_DstarDATA.pdf");
	if (dosp) l.SaveAs("DstarData/raw_yield_pt_logscale_DstarDATA.png");

	//evaluates the N systematics for each bin
	//for(int i=0;i<n_pt_bins;i++){
	// fit_syst_error_bin(datafile, pt_bins[i], pt_bins[i+1]);
	// }
}
//=====================================================================
void eta_analysis(RooWorkspace& ws, int n, TString etafile, TString datafile){

	TString dir_name = ".";
	TFile* f_wei = new TFile(dir_name + etafile, "recreate");

	RooAbsPdf* model = ws.pdf("model");
  	RooDataSet* data = (RooDataSet*) ws.data("data");
	#if particle == 0
	RooRealVar* Dseta = ws.var("Dseta");
  	RooRealVar Dsmass = *(ws.var("Dsmass"));
	#elif particle == 1
	RooRealVar* D0Kpieta = ws.var("D0Kpieta");
  	RooRealVar D0Kpimass = *(ws.var("D0Kpimass"));
	#endif
  	
	double mean_w = 0.;
	double mean_s = 0.;

	const int n_eta_bins = 10;
  	double eta_bins [n_eta_bins + 1] = {0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.1};

	double eta_mean[n_eta_bins];
	double eta_low[n_eta_bins];
	double eta_high[n_eta_bins];

	double yield[n_eta_bins];
	double yield_err_low[n_eta_bins];
	double yield_err_high[n_eta_bins];
	double yield_err_syst[n_eta_bins];
	double m_yield_err_syst[n_eta_bins];

	for(int k = 1; k<n_eta_bins; k++){ yield_err_syst[k] = 0; m_yield_err_syst[k] = 0; }

	RooDataSet* data_eta, data_w, data_wp;
	RooFitResult* fit_eta;
	RooRealVar* n_sig_eta;
	RooRealVar* n_comb_eta;

	const int n_pdf_syst=4;
	//number of pdfs
	TString syst_src[n_pdf_syst]={"nominal","bkg_poly","signal1gauss","bkg_range"};
	//array of pdfs to evaluate the systematic uncertainty of the N signal
	double yield_syst[n_eta_bins][n_pdf_syst]; 
	//value of the systematic uncertainty

	for(int i=0;i<n_eta_bins;i++)
	{	
		//select data subset corresponding to Eta bin
		#if particle == 0
		data_eta = (RooDataSet*) data->reduce(Form("fabs(Dseta)>%lf",eta_bins[i]));
		data_eta = (RooDataSet*) data_eta->reduce(Form("fabs(Dseta)<%lf",eta_bins[i+1]));
		#elif particle == 1
		data_eta = (RooDataSet*) data->reduce(Form("fabs(D0Kpieta)>%lf",eta_bins[i]));
		data_eta = (RooDataSet*) data_eta->reduce(Form("fabs(D0Kpieta)<%lf",eta_bins[i+1]));
		#endif
		
		ws.import(*data_eta, Rename(Form("data_eta_%d",i)));
   
		//perform fit and save result
		fit_eta = model->fitTo(*data_eta, Minos(true), Save());

		//plots the fit result
		#if particle == 0
		RooPlot* massframe = Dsmass.frame(Title(""));
		#elif particle == 1
		RooPlot* massframe = D0Kpimass.frame(Title(""));
		#endif

		data_eta->plotOn(massframe);
		model->paramOn(massframe,Layout(0.60,0.99,0.89));
		model->plotOn(massframe, Range("all"));
	
		TCanvas b;
		massframe->Draw();
		TLatex* tex11 = new TLatex(0.15,0.75, Form("%.2f<|#eta|<%.2f",eta_bins[i],eta_bins[i+1]));
		tex11->SetNDC(kTRUE);
		tex11->SetLineWidth(2);
		tex11->SetTextSize(0.04);
		tex11->Draw();
		tex11 = new TLatex(0.15,0.8," (pp) 13 TeV");
		tex11->SetNDC(kTRUE);
		tex11->SetTextFont(42);
		tex11->SetLineWidth(2);
		tex11->SetTextSize(0.04);
		tex11->Draw();
		tex11 = new TLatex(0.15,0.85,"Preliminary Study");
		tex11->SetNDC(kTRUE);
		tex11->SetTextFont(42);
		tex11->SetTextSize(0.04);
		tex11->SetLineWidth(2);
		tex11->Draw();

		int lowbins = eta_bins[i]*10;
		int highbins = eta_bins[i+1]*10;
#if particle == 0
	#if Dataset == 0
		b.SaveAs(Form("DstarData/ETA%i_%iDstarDATA.png", lowbins, highbins));
		b.SaveAs(Form("DstarData/ETA%i_%iDstarDATA.pdf", lowbins, highbins));
	#elif Dataset == 1
		b.SaveAs(Form("DstarMC/ETA%i_%iDstarMC.png", lowbins, highbins));
		b.SaveAs(Form("DstarMC/ETA%i_%iDstarMC.pdf", lowbins, highbins));
	#elif Dataset == 2
		b.SaveAs(Form("DstarMinBias/ETA%i_%iDstarMinBias.png", lowbins, highbins));
		b.SaveAs(Form("DstarMinBias/ETA%i_%iDstarMinBias.pdf", lowbins, highbins));
	#endif
#elif particle == 1
	#if Dataset == 0
		b.SaveAs(Form("D0Data/ETA%i_%iD0DATA.png", lowbins, highbins));
		b.SaveAs(Form("D0Data/ETA%i_%iD0DATA.pdf", lowbins, highbins));
	#elif Dataset == 1
		b.SaveAs(Form("D0MC/ETA%i_%iD0MC.png", lowbins, highbins));
		b.SaveAs(Form("D0MC/ETA%i_%iD0MC.pdf", lowbins, highbins));
	#elif Dataset == 2
		b.SaveAs(Form("D0MinBias/ETA%i_%iD0MinBias.png", lowbins, highbins));
		b.SaveAs(Form("D0MinBias/ETA%i_%iD0MinBias.pdf", lowbins, highbins));
	#endif		
#endif   
	}
}
//=====================================================================
void ForLifetime(RooWorkspace& ws){
	cout << "---------------------------------" << endl;
	cout << "ForLifetime" << endl;
	cout << "---------------------------------" << endl;

#if particle == 0
	RooRealVar Dslifetime = *(ws.var("Dslifetime"));
#elif particle == 1
	RooRealVar D0lifetime = *(ws.var("D0lifetime"));
#endif

	//Get the dataset with Wights
	RooDataSet* data = (RooDataSet*) ws.data("data");
	RooDataSet* dataW = (RooDataSet*) ws.data("dataWithSWeights");
	//Get the from the "dataWithSWeights" dataser only the signal
	RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
	RooDataSet* dataWBp2;

#if particle == 0
	dataWBp2 = (RooDataSet*) dataWBp->reduce("Dslifetime>(0.0*pow(10,-12))");
#elif particle == 1
	dataWBp2 = (RooDataSet*) dataWBp->reduce("D0lifetime>(0.0*pow(10,-12))");
#endif

	ws.import(*dataWBp2, Rename("data_lifetime"));
	//------------------------------------------------
	// P D F    L I F E T I M E
	//------------------------------------------------
	Double_t lambdaNonPrompt = -1./(1.520*pow(10,-12)); //B decay
	Double_t lambdaPrompt = -1./(4.101*pow(10,-13)); //D0 decay

	RooRealVar lambda1("lambda1", "lambda1", lambdaPrompt);
	RooRealVar lambda2("lambda2","lambda2", lambdaNonPrompt);

#if particle == 0
	RooExponential lifetime1("lifetime1", "lifetime1", Dslifetime, lambda1);
  	RooExponential lifetime2("lifetime2", "lifetime2", Dslifetime, lambda2);
#elif particle == 1
	RooExponential lifetime1("lifetime1", "lifetime1", D0lifetime, lambda1);
  	RooExponential lifetime2("lifetime2", "lifetime2", D0lifetime, lambda2);
#endif
  	
   lambda1.setConstant();
   lambda2.setConstant();
	RooRealVar prompt("prompt","prompt", 1000., 0.0, 5000000.);
	RooRealVar non_prompt("non_prompt","non_prompt", 1000., 0.0, 5000000.);
	RooAddPdf LifetimeModel("LifetimeModel","LifetimeModel", RooArgList(lifetime1,lifetime2), RooArgList(prompt,non_prompt));

	ws.import(LifetimeModel);

#if particle == 0
	#if Dataset == 0
	Dslifetime.setRange("range1",0.3*pow(10,-12),1.6*pow(10,-12));
	#elif Dataset == 1
	Dslifetime.setRange("range1",0.3*pow(10,-12),0.9*pow(10,-12));
	#elif Dataset == 2
	Dslifetime.setRange("range1",0.3*pow(10,-12),0.9*pow(10,-12));
	#endif
#elif particle == 1
	#if Dataset == 0
	D0lifetime.setRange("range1",0.8*pow(10,-12),5.0*pow(10,-12));
	#elif Dataset == 1
	D0lifetime.setRange("range1",0.6*pow(10,-12),2.8*pow(10,-12));
	#elif Dataset == 2
	D0lifetime.setRange("range1",0.7*pow(10,-12),5.0*pow(10,-12));
	#endif		
#endif

	//FIT the data with the model
	//LifetimeModel.fitTo(*dataWBp2,Range("all"));
	//LifetimeModel.fitTo(*dataWBp2,Extended());
	LifetimeModel.fitTo(*dataWBp2,Range("range1"));
	
	TCanvas d;
  	d.SetTitle("");
#if particle == 0
	RooPlot* massframe = Dslifetime.frame();
	massframe->SetTitle("D0(from D*) Lifetime");
#elif particle == 1
	RooPlot* massframe = D0lifetime.frame();
	massframe->SetTitle("D0 Lifetime");
#endif

	//d.SetLogy();
	//massframe->GetYaxis()->SetLogy();
	//massframe->SetLogy();

	dataWBp2->plotOn(massframe, RooFit::Name("dataWBp2") );
	LifetimeModel.plotOn(massframe, RooFit::Name("Fit"),Range("range1"),LineColor(kRed),LineStyle(1),LineWidth(2));
	LifetimeModel.plotOn(massframe, RooFit::Name("lifetime1"),Components("lifetime1"),Range("range1"),LineColor(kBlue),LineStyle(kDashed));
	LifetimeModel.plotOn(massframe, RooFit::Name("lifetime2"),Components("lifetime2"),Range("range1"),LineColor(kOrange),LineStyle(kDashed));
	LifetimeModel.paramOn(massframe,Layout(0.75,0.99,0.75));
	massframe->getAttText()->SetTextSize(0.028);
	massframe->GetYaxis()->SetTitleOffset(1.3);
	massframe->SetXTitle("Lifetime [s]");
	massframe->Draw();

	double chis = massframe->chiSquare();
	TLatex* tex13 = new TLatex(0.15, 0.85, Form("#chi^{2}/DOF = %.3lf",chis));
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04);	
	tex13->Draw();

	double totalEntries = dataWBp2->sumEntries();

	//TLegend *leg = new TLegend (0.4, 0.5, 0.6, 0.7); 
	TLegend *leg = new TLegend (0.45,0.70,0.65,0.85);
	#if Dataset == 0
	leg->AddEntry(massframe->findObject("dataWBp2"), "Data", "LP");
	#elif Dataset == 1
	leg->AddEntry(massframe->findObject("dataWBp2"), "MC", "LP");
	#elif Dataset == 2
	leg->AddEntry(massframe->findObject("dataWBp2"), "MinBias", "LP");
	#endif
	leg->AddEntry("dataWBp2",(Form("Entries: %2.0f", totalEntries))," ");
	leg->AddEntry("Fit","Sig + Bkg","L");
	leg->AddEntry("lifetime1","Prompt","L");
	leg->AddEntry("lifetime2","Non-Prompt","L");
	leg->Draw();

#if particle == 0
	#if Dataset == 0
	d.SaveAs("DstarData/DsLifetimeDATA.pdf");
	#elif Dataset == 1
	d.SaveAs("DstarMC/DsLifetimeMC.pdf");
	#elif Dataset == 2
	d.SaveAs("DstarMinBias/DsLifetimeMinBias.pdf");
	#endif
#elif particle == 1
	#if Dataset == 0
	d.SaveAs("D0DATA/D0LifetimeDATA.pdf");
	#elif Dataset == 1
	d.SaveAs("D0MC/D0LifetimeMC.pdf");
	#elif Dataset == 2
	d.SaveAs("D0MinBias/D0LifetimeMinBias.pdf");
	#endif
#endif

	cout << "---------------------------------" << endl;
	cout << "End ForLifetime" << endl;
	cout << "---------------------------------" << endl;	
}

