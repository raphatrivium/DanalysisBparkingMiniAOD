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
#include "RooGaussModel.h"
#include "RooDecay.h"

using namespace RooStats;
using namespace RooFit;
using namespace std;

// see below for implementation

void set_up_workspace_variables(RooWorkspace& ws);
void AddData(RooWorkspace& ws, TString f_input);
void build_pdf(RooWorkspace& w, std::string choice = "nominal");
void AddModel(RooWorkspace& ws);
void plot_complete_fit(RooWorkspace& ws);
std::vector<TH1D*> sideband_subtraction(RooWorkspace* w, int* n, int n_var);
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n, bool sPLOTlMethod); 
void DoSPlot(RooWorkspace& ws);
void MakePlots(RooWorkspace& ws, int n, TString label);
void pT_analysis(RooWorkspace& ws, int n, TString ptfile, TString datafile);
void eta_analysis(RooWorkspace& ws, int n, TString ptfile, TString datafile);
void ForLifetime(RooWorkspace& ws);
void GenLifetime(RooWorkspace& ws);
void ForLifetimeV2(RooWorkspace& ws);

//--------------------------------------
// D E F I N E  T H E  P A R T I C L E
#define particle 0 // 0= D* ; 1= D0
//--------------------------------------
// D E F I N E  T H E  D A T A S E T
#define Dataset 0 // 0= Data ; 1= MC ; 2= MinBias
//--------------------------------------
// D E F I N E   M E T H O D   I N   L I F E T I M E   F I T
#define LifeSideBand 0 // 0= sPLOT ; 1= SideBand

void DstarD0Fit()
{	//For the out put of time taken to this program run
	clock_t tStart = clock();

	//-------------------------------------------------------------
	//D A T A S E Ts   O P T I O S
	//-------------------------------------------------------------
	#if Dataset == 0 // Data
	TString input_file_data = "BparkingDataset_Data.root"; // BparkingDataset_Data.root BparkingDatasetSig1.root
	#elif Dataset == 1 // MC
	TString input_file_data = "MC_DStarToD0Pi_D0KPi.root";
	#elif Dataset == 2 // MinBias
	TString input_file_data = "MC_MinBias.root";
	#endif
	//-------------------------------------------------------------
	//V A R I A B L E S   A N D   B I N S   F O R   E A C H   P A R T I C L E
	//-------------------------------------------------------------
	//variables are the branches in input RootFile
#if particle == 0 // Dstar
	TString variables[] = {"Dseta","Dsphi","Dspt","D0fromDSs3D",
									"D0mass","D0eta","D0phi","D0pt", "Dslifetime", "DstarCt"};
	//Number of varibles in  variables[]
	int n_var = sizeof(variables)/sizeof(variables[0]);
	//Number of bins
	int n_bins[]= {100, 100, 10, 10, 20, 10, 10, 10, 10, 10};
#elif particle == 1 // D0
	TString variables[] = {"D0Kpipt","D0Kpieta","D0Kpiphi","D0lifetime", "D0Ct"};
	int n_var = sizeof(variables)/sizeof(variables[0]);
	int n_bins[]= {100, 100, 10, 10, 20, 10, 10, 10, 10, 10, 15};
#endif
	//-------------------------------------------------------------
	//F U N C T I O N S
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
	//build_pdf(*MyWS);
	//-------------------------------------------------------------
	//Fit Signal + Background. (depends of AddData())
	plot_complete_fit(*MyWS);
	//-------------------------------------------------------------	
	//sideband_sub histograms - Returns a vector of TH1D
	//std::vector<TH1D*> histos_sideband_sub;
	//histos_sideband_sub = sideband_subtraction(MyWS, n_bins, n_var);
	//-------------------------------------------------------------
	// do sPlot.This wil make a new dataset with sWeights added for every event. (depends of plot_complete_fit())
	DoSPlot(*MyWS);
	//-------------------------------------------------------------
	// Make some plots showing the discriminating variable. (depends of DoSPlot())
	//for (int i=0; i < n_var; i++ ) {MakePlots(*MyWS, n_var, variables[i]);}
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
	ForLifetimeV2(*MyWS); //Correctly way
	#if Dataset == 0 // Data
	#elif Dataset == 1 // MC
	//GenLifetime(*MyWS); //(depends of AddData())
	#elif Dataset == 2 // MinBias
	#endif	
	//-------------------------------------------------------------
  	//cleanup
  	delete MyWS;
	//-------------------------------------------------------------
	//T I M E   O U T P U T
	//-------------------------------------------------------------
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	printf("          : %.2fm\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC*60.));
	printf("          : %.2fh\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC*60.*60.));		
	cout << "---------------------------------" << endl;
	cout << "E N D   P R O G R A M" << endl;
	cout << "---------------------------------" << endl;
}//end main program
//=====================================================================
//=====================================================================
void set_up_workspace_variables(RooWorkspace& ws)
{
	cout << "---------------------------------" << endl;
	cout << "set_up_workspace_variables" << endl;
	cout << "---------------------------------" << endl;
#if particle == 0 // Dstar
	//------------------------------------------
	//set range of observable
	//------------------------------------------
	Double_t lowRange = 1.93, highRange = 2.10;
	Double_t lowRange2 = 1.77, highRange2 = 1.95;
	Double_t lowEta = -4., highEta = 4.;
	Double_t lowPhi = -4., highPhi = 4.;
	Double_t lowdxy = -0.1, highdxy = 0.1;
	Double_t lowdz = -0.1, highdz = 0.1;
	Double_t lowlifetime = 0.*(pow(10,-12)); 	Double_t highlifetime = 1.6*(pow(10,-12));
	Double_t lowProperTime = 0.0; 	Double_t highProperTime = 0.07;
	//------------------------------------------
	//make a RooRealVar for the observables
	//------------------------------------------
	RooRealVar Dsmass("Dsmass","D* M_{inv}",lowRange,highRange, "GeV/c^{2}");
	ws.import(Dsmass); // import to your worksapce
	RooRealVar Dseta("Dseta","#eta",lowEta,highEta); ws.import(Dseta);
	RooRealVar Dsphi("Dsphi","#phi",lowPhi,highPhi); ws.import(Dsphi);
	RooRealVar Dspt("Dspt","Transverse Momentum",0.,100.,"GeV/c");	ws.import(Dspt);
	RooRealVar D0fromDSs3D("D0fromDSs3D","D0fromDSs3D",0.,5.); ws.import(D0fromDSs3D);
	RooRealVar D0mass("D0mass","D0mass",lowRange2,highRange2); ws.import(D0mass);
	RooRealVar D0phi("D0phi","D0phi",lowPhi,highPhi); ws.import(D0phi);
	RooRealVar D0eta("D0eta","D0eta",lowEta,highEta); ws.import(D0eta);
	RooRealVar D0pt("D0pt","D0pt",4.,100.); ws.import(D0pt); 
	RooRealVar Dslifetime("Dslifetime","Dslifetime",lowlifetime,highlifetime); ws.import(Dslifetime);
	RooRealVar DstarCt("DstarCt","DstarCt",lowProperTime,highProperTime); ws.import(DstarCt);
	#if Dataset == 0 // Data
	#elif Dataset == 1 // MC
	RooRealVar PUWeight("PUWeight","PUWeight",0.,10.,"Weight"); ws.import(PUWeight);
	//RooRealVar MCD0lifetime("MCD0lifetime","MCD0lifetime",lowlifetime,highlifetime);
	//ws.import(MCD0lifetime);
	RooRealVar MCD0lifetimeMatching("MCD0lifetimeMatching","MCD0lifetimeMatching",lowlifetime,highlifetime);	ws.import(MCD0lifetimeMatching);
	#elif Dataset == 2 // MinBias
	RooRealVar PUWeight("PUWeight","PUWeight",0.,10.,"Weight");	ws.import(PUWeight);
	// Contamination
	RooRealVar DstarFromBmass("DstarFromBmass","D* M_{inv}", lowRange, highRange, "GeV/c^{2}"); 	ws.import(DstarFromBmass);
	RooRealVar DstarFromBpt("DstarFromBpt","Transverse Momentum",0.,100.,"GeV/c"); ws.import(DstarFromBpt);
	RooRealVar DstarFromBeta("DstarFromBeta","#eta", lowEta, highEta); ws.import(DstarFromBeta);
	RooRealVar DstarFromPPmass("DstarFromPPmass","D* M_{inv}", lowRange, highRange, "GeV/c^{2}"); ws.import(DstarFromPPmass);
	RooRealVar DstarFromPPpt("DstarFromPPpt","Transverse Momentum",0.,100.,"GeV/c");	ws.import(DstarFromPPpt);
	RooRealVar DstarFromPPeta("DstarFromPPeta","#eta", lowEta, highEta);	ws.import(DstarFromPPeta);
	#endif
#elif particle == 1 // D0
  	// set range of observable
	Double_t lowRange = 1.77, highRange = 1.96;
	Double_t lowEta = -4., highEta = 4.;
	Double_t lowPhi = -4., highPhi = 4.;
	Double_t lowdxy = -0.1, highdxy = 0.1;
	Double_t lowdz = -0.1, highdz = 0.1;
	Double_t lowlifetime = 0.0*(pow(10,-12));
	Double_t highlifetime = 5.0*(pow(10,-12));
	// make a RooRealVar for the observables
	RooRealVar D0Kpimass("D0Kpimass","D0 Invariant Mass",lowRange,highRange, "GeV/c^{2}");	ws.import(D0Kpimass);
	RooRealVar D0Kpieta("D0Kpieta","#eta",lowEta,highEta); ws.import(D0Kpieta);
	RooRealVar D0Kpiphi("D0Kpiphi","#phi",lowPhi,highPhi); ws.import(D0Kpiphi);
	RooRealVar D0Kpipt("D0Kpipt","Transverse Momentum",0.,100.,"GeV/c");	ws.import(D0Kpipt);
	RooRealVar D0lifetime("D0lifetime","D0lifetime",lowlifetime,highlifetime,"s"); ws.import(D0lifetime);
	#if Dataset == 0 // Data
	#elif Dataset == 1 // MC
	RooRealVar PUWeight("PUWeight","PUWeight",0.,20.,"Weight"); ws.import(PUWeight);
	//RooRealVar MCpromptD0lifetime("MCpromptD0lifetime","MCpromptD0lifetime",lowlifetime,highlifetime); ws.import(MCpromptD0lifetime);
	RooRealVar MCpromptD0lifetimeMatching("MCpromptD0lifetimeMatching","MCpromptD0lifetimeMatching",lowlifetime,highlifetime);	ws.import(MCpromptD0lifetimeMatching);
	#elif Dataset == 2 // MinBias
	RooRealVar PUWeight("PUWeight","PUWeight",0.,20.,"Weight"); ws.import(PUWeight);
	// Contamination
	RooRealVar D0FromBmass("D0FromBmass","D* M_{inv}", lowRange, highRange, "GeV/c^{2}");	ws.import(D0FromBmass);
	RooRealVar D0FromBpt("D0FromBpt","Transverse Momentum",0.,100.,"GeV/c"); ws.import(D0FromBpt);
	RooRealVar D0FromBeta("D0FromBeta","#eta", lowEta, highEta); ws.import(D0FromBeta);
	RooRealVar D0FromPPmass("D0FromPPmass","D* M_{inv}", lowRange, highRange, "GeV/c^{2}"); ws.import(D0FromPPmass);
	RooRealVar D0FromPPpt("D0FromPPpt","Transverse Momentum",0.,100.,"GeV/c");	ws.import(D0FromPPpt);
	RooRealVar D0FromPPeta("D0FromPPeta","#eta", lowEta, highEta);	ws.import(D0FromPPeta);
	#endif
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
#if particle == 0 // Dstar
	TTree* t1_data = (TTree *)file->Get("t_analysis");  //   t_D0MCMatching
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
	arg_list.add(*(ws.var("DstarCt")));
	#if Dataset == 0 // Data
	RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);
	ws.import(*data, Rename("data"));
	#elif Dataset == 1 // MC
	//define weights for MC
	arg_list.add(*(ws.var("PUWeight")));
	RooDataSet* data = new RooDataSet("data", "data", t1_data, arg_list);
	RooRealVar PUWeight = *(ws.var("PUWeight"));
	// Construct formula to calculate weight for events
	RooFormulaVar wFunc("w","event weight","1.*PUWeight",PUWeight) ;
	// Add column with variable w to previously generated dataset
	//RooRealVar* w = (RooRealVar*) data->addColumn(PUWeight);
	RooRealVar* w = (RooRealVar*) data->addColumn(wFunc);
	// Instruct dataset wdata in interpret w as event weight rather than as observable
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName());
	wdata.Print();
	data->get(0)->Print("Dsmass"); // For printing fit variables
	data->get(1)->Print("Dsmass");
	ws.import(wdata, Rename("data")); //*data  wdata

	//For lifetime of the generated particle
	TTree* t2_data = (TTree *)file->Get("t_DsMatching"); // t_DsMC    t_DsMatching
	RooArgList arg_list2 ("arg_list2");
	//arg_list2.add(*(ws.var("MCD0lifetime")));
	arg_list2.add(*(ws.var("MCD0lifetimeMatching")));
	RooDataSet* data2 = new RooDataSet("data2", "data2", t2_data, arg_list2);
	ws.import(*data2, Rename("data2"));

	#elif Dataset == 2 // MinBias
	arg_list.add(*(ws.var("PUWeight")));
	RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);
	RooRealVar PUWeight = *(ws.var("PUWeight"));
	RooFormulaVar wFunc("w","event weight","1.*PUWeight",PUWeight) ;
	RooRealVar* w = (RooRealVar*) data->addColumn(wFunc);
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName());
	ws.import(wdata, Rename("data"));
	//contamination
	TTree* t_Dstarcont = (TTree *)file->Get("t_DstarContamination");
	RooArgList cont_list1 ("cont_list1");
	cont_list1.add(*(ws.var("DstarFromBmass")));
	cont_list1.add(*(ws.var("DstarFromBpt")));
	cont_list1.add(*(ws.var("DstarFromBeta")));
	RooDataSet* data_contB = new RooDataSet("data_contB", "data_contB", t_Dstarcont, cont_list1);
	ws.import(*data_contB, Rename("data_contB"));
	//----------------------------------------------------
	RooArgList cont_list2 ("cont_list2");
	cont_list2.add(*(ws.var("DstarFromPPmass")));
	cont_list2.add(*(ws.var("DstarFromPPpt")));
	cont_list2.add(*(ws.var("DstarFromPPeta")));
	RooDataSet* data_contPP = new RooDataSet("data_contPP", "data_contPP", t_Dstarcont, cont_list2);
	ws.import(*data_contPP, Rename("data_contPP"));
	#endif

#elif particle == 1 // D0
	TTree* t1_data = (TTree *)file->Get("t_D0analysis");
	RooArgList arg_list ("arg_list");

	arg_list.add(*(ws.var("D0Kpimass")));
	arg_list.add(*(ws.var("D0Kpieta")));
	arg_list.add(*(ws.var("D0Kpiphi")));
	arg_list.add(*(ws.var("D0Kpipt")));
	arg_list.add(*(ws.var("D0lifetime")));
	arg_list.add(*(ws.var("D0Ct")));
	#if Dataset == 0 // Data
	#elif Dataset == 1 // MC
	RooRealVar PUWeight = *(ws.var("PUWeight"));
	#elif Dataset == 2 // MinBias
	RooRealVar PUWeight = *(ws.var("PUWeight"));
	#endif

	RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);
	#if Dataset == 0 // Data
	ws.import(*data, Rename("data"));
	#elif Dataset == 1 // MC
	RooRealVar* w = (RooRealVar*) data->addColumn(PUWeight) ;
	RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName()) ;
	ws.import(wdata, Rename("data"));

	//For lifetime of the generated particle
	TTree* t2_data = (TTree *)file->Get("t_D0Matching"); //t_D0MC   t_D0Matching   
	RooArgList arg_list2 ("arg_list2");
	//arg_list2.add(*(ws.var("MCpromptD0lifetime")));
	arg_list2.add(*(ws.var("MCpromptD0lifetimeMatching")));
	RooDataSet* data2 = new RooDataSet("data2", "data2", t2_data, arg_list2);
	ws.import(*data2, Rename("data2"));
	
	#elif Dataset == 2 // MinBias
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
//=====================================================================
void AddModel(RooWorkspace& ws) {
	cout << "---------------------------------" << endl;
	cout << "AddModel" << endl;
	cout << "---------------------------------" << endl;

#if particle == 0 // Dstar
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

	#if Dataset == 0 // Data
	sigma1.setConstant();
	sigma2.setConstant();
	#elif Dataset == 1 // MC
	#elif Dataset == 2 // MinBias
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

#elif particle == 1 // D0
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
//=====================================================================
void plot_complete_fit(RooWorkspace& ws){
	cout << "---------------------------------" << endl;
	cout << "plot_complete_fit" << endl;
	cout << "---------------------------------" << endl;

	RooAbsPdf* model = ws.pdf("model");
	//double fwhmModel= model->get_fwhm();
	RooDataSet* data = (RooDataSet*) ws.data("data");
	RooRealVar* sigma1 = ws.var("sigma1");
	RooRealVar* sigma2 = ws.var("sigma2");
#if particle == 0 // Dstar
	RooRealVar Dsmass = *(ws.var("Dsmass"));
	RooPlot* massframe = Dsmass.frame();
	massframe->SetTitle("D* Invariant Mass");
#elif particle == 1 // D0
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
	data->plotOn(massframe, RooFit::Name("Data") );
	model->plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
	//OBS!!! if you get the "chiSquare()" after plot other components its value may change
	double chis = massframe->chiSquare();
	model->plotOn(massframe, RooFit::Name("Signal"),Components("sig"),Range("all"), LineColor(kOrange), LineStyle(kDashed), DrawOption("F"), FillColor(kOrange), FillStyle(2015));
	model->plotOn(massframe, RooFit::Name("gauss1"),Components("gauss1"),Range("all"), LineColor(kGreen), LineStyle(kDashed));
	model->plotOn(massframe, RooFit::Name("gauss2"),Components("gauss2"),Range("all"), LineColor(kGreen), LineStyle(kDashed));
	model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("bkg"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
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
	//-------------------------------
	//L E G E N D
	//----------------------------------
	#if Dataset == 0 // Data
	TLatex* tex11 = new TLatex(0.6,0.8,"Data (pp) 13 TeV");
	#elif Dataset == 1 // MC
	TLatex* tex11 = new TLatex(0.6,0.8,"MC (pp) 13 TeV");
	#elif Dataset == 2 // MinBias
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

	

	/*TLatex* tex12 = new TLatex(0.15, 0.85, Form("#sigma_{exp} = %.3lf #pm %.3lf",sigma1_str,sigma1_err));
	tex12->SetNDC(kTRUE); tex12->SetTextFont(42); tex12->SetTextSize(0.04);
	tex12->Draw();*/

	TLatex* tex13 = new TLatex(0.15, 0.85, Form("#chi^{2}/DOF = %.3lf",chis));
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04);
	tex13->Draw();

	double totalEntries = data->sumEntries();

#if particle == 0 // Dstar
	#if Dataset == 0 // Data
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "Data", "LP");
	#elif Dataset == 1 // MC
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "MC", "LP");
	#elif Dataset == 2 // MinBias
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "MinBias", "LP");
	#endif
#elif particle == 1 // D0
	#if Dataset == 0 // Data
	TLegend *leg = new TLegend (0.15,0.25,0.35,0.45);
	leg->AddEntry(massframe->findObject("Data"), "Data", "LP");
	#elif Dataset == 1 // MC
	TLegend *leg = new TLegend (0.15,0.55,0.35,0.75);
	leg->AddEntry(massframe->findObject("Data"), "MC", "LP");
	#elif Dataset == 2 // MinBias
	TLegend *leg = new TLegend (0.15,0.25,0.35,0.45);
	leg->AddEntry(massframe->findObject("Data"), "MinBias", "LP");
	#endif		
#endif

	leg->AddEntry("Data",(Form("Entries: %2.0f", totalEntries))," ");
	leg->AddEntry("Fit","Sig + Bkg","L");
	leg->AddEntry("Signal","Signal","LP");
	leg->AddEntry("Combinatorial","Combinatorial","LP");
	leg->Draw();

	//Pull dists
	RooHist* pull_hist = massframe->pullHist("Data","Fit");
#if particle == 0 // Dstar
	RooPlot *pull_plot = Dsmass.frame();
#elif particle == 1 // D0
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

#if particle == 0 // Dstar
	#if Dataset == 0 // Data
		d.SaveAs("DstarData/FitDstarDATA.pdf");
		d.SaveAs("DstarData/FitDstarDATA.png");
	#elif Dataset == 1 // MC
		d.SaveAs("DstarMC/FitDstarMC.pdf");
		d.SaveAs("DstarMC/FitDstarMC.png");
	#elif Dataset == 2 // MinBias
		d.SaveAs("DstarMinBias/FitDstarMinBias.pdf");
		d.SaveAs("DstarMinBias/FitDstarMinBias.png");
	#endif
#elif particle == 1 // D0
	#if Dataset == 0 // Data
		d.SaveAs("D0DATA/FitD0DATA.pdf");
		d.SaveAs("D0DATA/FitD0DATA.png");
	#elif Dataset == 1 // MC
		d.SaveAs("D0MC/FitD0MC.pdf");
		d.SaveAs("D0MC/FitD0MC.png");
	#elif Dataset == 2 // MinBias
		d.SaveAs("D0MinBias/FitD0MinBias.pdf");
		d.SaveAs("D0MinBias/FitD0MinBias.png");
	#endif		
#endif

	cout << "---------------------------------" << endl;
	cout << "End plot_complete_fit" << endl;
	cout << "---------------------------------" << endl;
}//End plot_complete_fit
//=====================================================================
//=====================================================================
std::vector<TH1D*> sideband_subtraction(RooWorkspace* ws, int* nBins, int n_var){
	cout << "---------------------------------" << endl;
	cout << "sideband_subtraction()" << endl;
	cout << "---------------------------------" << endl;
	RooDataSet* data = (RooDataSet*) ws->data("data");
	RooAbsPdf* fit_side = ws->pdf("bkg");
	RooDataSet* data_side;
	RooDataSet* data_central;
	
	vector<RooRealVar> variables;								
	

#if particle == 0 // Dstar
	variables.push_back(*(ws->var("Dsmass")));
	variables.push_back(*(ws->var("Dspt")));
	variables.push_back(*(ws->var("Dslifetime")));
	variables.push_back(*(ws->var("Dseta")));
	variables.push_back(*(ws->var("Dsphi")));
	variables.push_back(*(ws->var("D0fromDSs3D")));
	variables.push_back(*(ws->var("D0mass")));
	variables.push_back(*(ws->var("D0eta")));
	variables.push_back(*(ws->var("D0phi")));
	variables.push_back(*(ws->var("D0pt")));
#elif particle == 1 // D0
	variables.push_back(*(ws->var("D0Kpimass")));
	variables.push_back(*(ws->var("D0Kpipt")));
	variables.push_back(*(ws->var("D0lifetime")));
	variables.push_back(*(ws->var("D0Kpieta")));
	variables.push_back(*(ws->var("D0Kpiphi")));
#endif

	//-------------------------------------------------------
	//D E F I N I N G   D A T A S E T S
	//-------------------------------------------------------
#if particle == 0 // Dstar
	//double left = particle ? 1.98 : 1.77; 
	//double right = particle ? 2.05 : 1.96; 
	double left = 1.99; //1.99    2.01032-2.355*0.0189
	double right = 2.035; //2.035   2.01032+2.355*0.0189
	//data_side = particle ? (RooDataSet*)data->reduce(Form("Dsmass>%lf or Dsmass<%lf", left, right)) : (RooDataSet*)data->reduce(Form("Dsmass>%lf",right));
	//BACKGROUND
	data_side = (RooDataSet*)data->reduce(Form("Dsmass>%lf || Dsmass<%lf", right, left));
	//SIGNAL
	data_central = (RooDataSet*)data->reduce(Form("Dsmass>%lf",left));
	data_central = (RooDataSet*)data_central->reduce(Form("Dsmass<%lf",right));
#elif particle == 1 // D0
	double left = 1.84; //1.99    2.01032-2.355*0.0189
	double right = 1.89; //2.035   2.01032+2.355*0.0189
	//BACKGROUND
	data_side = (RooDataSet*)data->reduce(Form("D0Kpimass>%lf || D0Kpimass<%lf", right, left));
	//SIGNAL
	data_central = (RooDataSet*)data->reduce(Form("D0Kpimass>%lf",left));
	data_central = (RooDataSet*)data_central->reduce(Form("D0Kpimass<%lf",right));
#endif
	
	//-------------------------------------------------------
	//INTEGRATING THE BACKGROUND DISTRIBUTION
	//-------------------------------------------------------
	RooAbsReal* int_right = fit_side->createIntegral(variables[0], "right");
	RooAbsReal* int_left = fit_side->createIntegral(variables[0], "left");
	RooAbsReal* int_peak = fit_side->createIntegral(variables[0], "peak");
	cout << "Integral right band: " << int_right->getVal() << endl;
	cout << "Integral left band : " << int_left->getVal() << endl;
	//double factor = particle ? (int_peak->getVal())/(int_right->getVal() + int_left->getVal()) : (int_peak->getVal())/(int_right->getVal());
	//double factor = (int_peak->getVal()) / (int_right->getVal() + int_left->getVal());
	double factor = 1.;
	cout << "Factor: " << factor << endl;

	for(int i=0; i<n_var; i++){std::cout << "bins: " << nBins[i] << std::endl;}

	std::vector<TH1D*> histos;
#if particle == 0 // Dstar
	histos.push_back(create_histogram(variables[0],"SideBand_Dsmass", factor, data_side, data_central, data, nBins[0], false));
	histos.push_back(create_histogram(variables[1],"SideBand_Dspt", factor, data_side, data_central, data, nBins[0], false) );
	histos.push_back(create_histogram(variables[2],"SideBand_Dslifetime", factor, data_side, data_central, data, nBins[1], false));
#elif particle == 1 // D0
	histos.push_back(create_histogram(variables[0],"SideBand_D0Kpimass", factor, data_side, data_central, data, nBins[0], false));
	histos.push_back(create_histogram(variables[1],"SideBand_D0Kpipt", factor, data_side, data_central, data, nBins[0], false) );
	histos.push_back(create_histogram(variables[2],"SideBand_D0lifetime", factor, data_side, data_central, data, nBins[1], false));
#endif
	
	return histos;

	ws->import(*data_central, Rename("data_central"));
	
	cout << "---------------------------------" << endl;
	cout << "END sideband_subtraction" << endl;
	cout << "---------------------------------" << endl;
}
//=====================================================================
//=====================================================================
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n, bool sPLOTlMethod){
	cout << "---------------------------------" << endl;
	cout << "BEGIN create_histogram()" << endl;
	cout << "---------------------------------" << endl;
	std::cout<< "Number of Bins= "<< n <<std::endl;
	//------------------------------------
	// B A C K G R O U N D   H I S T O G R A M 
	//------------------------------------
	TH1D* hist_side = (TH1D*)reduced->createHistogram("hist_side",var, Binning(n, var.getMin(), var.getMax()));
	hist_side->SetMarkerColor(kRed);
	hist_side->SetLineColor(kRed);
	hist_side->SetNameTitle("hist_side", "");
	hist_side->Scale(factor);
	hist_side->SetStats(0); //without Statistic Box  
	//------------------------------------
	// S I G N A L   H I S T O G R A M
	//------------------------------------
	TH1D* hist_peak = (TH1D*)central->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));
	hist_peak->SetMarkerColor(kBlue);
	hist_peak->SetLineColor(kBlue);
	hist_peak->SetNameTitle("hist_peak", "");
	hist_peak->SetStats(0);
	//------------------------------------
	// S I G N A L + B A C K G R O U N D    H I S T O G R A M
	//------------------------------------
	//TH1D* hist_total = (TH1D*)total->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));
	TH1D* hist_total = new TH1D(*hist_peak);
	hist_total->Add(hist_side, factor);
	hist_total->SetMarkerColor(kBlack);
	hist_total->SetLineColor(kBlack);
	hist_total->SetNameTitle(var.GetName(), "");
	hist_total->SetStats(0);
	//------------------------------------
	// C R E A T I N G   C A N V A S
	//------------------------------------
	TCanvas c;
	hist_total->Draw();
	hist_side->Draw("same");
	hist_peak->Draw("same");

	hist_total->SetXTitle(var.GetName());
	hist_side->SetXTitle(var.GetName());
	hist_peak->SetXTitle(var.GetName());
	hist_peak->GetYaxis()->SetRangeUser(0, 1.3*hist_peak->GetMaximum());
	//------------------------------------------------------------
	//LEGEND
	//------------------------------------------------------------
	//TLatex* tex = new TLatex(0.6,0.8,"1.5 nb^{-1} (PbPb) 5.02 TeV");
	TLatex* tex = new TLatex(0.68,0.8,"(pp) 13 TeV");
	tex->SetNDC(kTRUE);
	tex->SetLineWidth(2);
	tex->SetTextSize(0.04);
	tex->Draw();
	if ( sPLOTlMethod ){tex = new TLatex(0.60,0.85,"sPLOT-CMS Preliminary");}
	else {tex = new TLatex(0.60,0.85,"SideBand-CMS Preliminary");}
	tex->SetNDC(kTRUE);
	tex->SetTextFont(42);
	tex->SetTextSize(0.04);
	tex->SetLineWidth(2);
	tex->Draw();

	TLegend *leg = new TLegend (0.7, 0.5, 0.85, 0.65);
	#if Dataset == 0 // Data
	leg->AddEntry("hist_total", "Data", "l*");
	#elif Dataset == 1 // MC
	leg->AddEntry("hist_total", "MC", "l*");
	#elif Dataset == 2 // MinBias
	leg->AddEntry("hist_total", "MinBias", "l*");
	#endif
	leg->AddEntry("hist_side", "Background", "l");
	leg->AddEntry("hist_peak", "Signal", "l");
	leg->Draw("same");

	std::cout<<"name: "<<var.GetName()<<std::endl;
	std::cout<<"histo name: "<<hist_total->GetName()<<std::endl;
	//------------------------------------------------------------
	//S A V E   A S   I M A G E
	//------------------------------------------------------------
#if particle == 0 // Dstar
	#if Dataset == 0 // Data
		c.SaveAs("DstarData/"+name+"DATA.pdf");
		c.SaveAs("DstarData/"+name+"DATA.png");
	#elif Dataset == 1 // MC
		c.SaveAs("DstarMC/"+name+"MC.pdf");
		c.SaveAs("DstarMC/"+name+"MC.png");
	#elif Dataset == 2 // MinBias
		c.SaveAs("DstarMinBias/"+name+"MinBias.pdf");
		c.SaveAs("DstarMinBias/"+name+"MinBias.png");
	#endif
#elif particle == 1 // D0
	#if Dataset == 0 // Data
		c.SaveAs("D0DATA/"+name+"DATA.pdf");
		c.SaveAs("D0DATA/"+name+"DATA.png");
	#elif Dataset == 1 // MC
		c.SaveAs("D0MC/"+name+"MC.pdf");
		c.SaveAs("D0MC/"+name+"MC.png");
	#elif Dataset == 2 // MinBias
		c.SaveAs("D0MinBias/"+name+"MinBias.pdf");
		c.SaveAs("D0MinBias/"+name+"MinBias.png");
	#endif		
#endif

	return hist_total;
	cout << "---------------------------------" << endl;
	cout << "END create_histogram()" << endl;
	cout << "---------------------------------" << endl;
}
//=====================================================================
//=====================================================================
void DoSPlot(RooWorkspace& ws){
	cout << "---------------------------------" << endl;
	cout << "Begin DoSPlot()" << endl;
	cout << "---------------------------------" << endl;

	//we need the fit and the dataset previously saved in the woorkspace
	RooDataSet* data = (RooDataSet*) ws.data("data");
   RooAbsPdf* model = ws.pdf("model");
	
	RooDataSet* data_side;
	RooDataSet* data_central;

#if particle == 0 // Dstar
	double left = 1.99; //1.99    2.01032-2.355*0.0189
	double right = 2.035; //2.035   2.01032+2.355*0.0189
	//SIGNAL
	data_central = (RooDataSet*)data->reduce(Form("Dsmass>%lf",left));
	data_central = (RooDataSet*)data_central->reduce(Form("Dsmass<%lf",right));
	//BACKGROUND
	data_side = (RooDataSet*)data->reduce(Form("Dsmass>%lf || Dsmass<%lf", right, left));
	ws.import(*data_central, Rename("data_central"));
#elif particle == 1 // D0
	double left = 1.99; //1.99    2.01032-2.355*0.0189
	double right = 2.035; //2.035   2.01032+2.355*0.0189
	//SIGNAL
	data_central = (RooDataSet*)data->reduce(Form("D0Kpimass>%lf",left));
	data_central = (RooDataSet*)data_central->reduce(Form("D0Kpimass<%lf",right));
	ws.import(*data_central, Rename("data_central"));
	//BACKGROUND
	data_side = (RooDataSet*)data->reduce(Form("D0Kpimass>%lf || D0Kpimass<%lf", right, left));
#endif
	
	
   
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
	
	//?
	RooMsgService::instance().setSilentMode(true);

	//add sWeights to dataset based on model and yield variables
	//sPlot class adds a new variable that has the name of the corresponding yield + "_sw".
   SPlot* sData = new SPlot("sData","An sPlot",*data, model, RooArgList(*D0Yield,*BgYield));

	cout << endl <<"Yield of Ds is "<< D0Yield->getVal() << ". From sWeights it is " << sData->GetYieldFromSWeight("n_signal") << endl;
	cout << "Yield of background is "<< BgYield->getVal() << ". From sWeights it is " << sData->GetYieldFromSWeight("n_combinatorial") <<  endl;

	// The first ten entries	
	for(Int_t i=0; i < 10; i++)
	{	std::cout << " Ds Weight " << sData->GetSWeight(i,"n_signal")
                << " Bg Weight " << sData->GetSWeight(i,"n_combinatorial")
                << " Total Weight " << sData->GetSumOfEventSWeight(i)
                << std::endl;}

	//the reweighted data is saved in the woorkspace
	ws.import(*data, Rename("dataWithSWeights"));

	cout << "---------------------------------" << endl;
	cout << "End DoSPlot" << endl;
	cout << "---------------------------------" << endl;
}//End DoSPlot
//=====================================================================
void MakePlots(RooWorkspace& ws, int n, TString label){
	cout << "---------------------------------" << endl;
	cout << "Begin MakePlots()" << endl;
	cout << "---------------------------------" << endl;
	// make our canvas
	TCanvas* cdata = new TCanvas("sPlot","sPlot", 800, 600);
	cdata->Divide(2,2);

	RooAbsPdf* model = ws.pdf("model");
	RooAbsPdf* sig = ws.pdf("sig");
	RooAbsPdf* bkg = ws.pdf("bkg");
	RooRealVar* variable = ws.var(label);

	#if particle == 0 // Dstar
	RooRealVar* Dsmass  = ws.var("Dsmass");
	#elif particle == 1 // D0
	RooRealVar* D0Kpimass  = ws.var("D0Kpimass");
	#endif	

	RooRealVar* BpYield = ws.var("n_signal");
	RooRealVar* BgYield = ws.var("n_combinatorial");

	double sigYield = BpYield->getVal();
	double bkgYield = BgYield->getVal();
	double nbin = 100;

	RooDataSet* data = (RooDataSet*) ws.data("data");

	cdata->cd(1);
#if particle == 0 // Dstar
	RooPlot* mframe = Dsmass->frame();
	mframe->SetTitle("D* Mass Distribuition");
#elif particle == 1 // D0
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
	RooPlot* frames = variable->frame();
	data->plotOn(frames);
	#if Dataset == 0 // Data
	frames->SetTitle(label + ": Total Sample (Data)");
	#elif Dataset == 1 // MC
	frames->SetTitle(label + ": Total Sample (MC)");
	#elif Dataset == 2 // MinBias
	frames->SetTitle(label + ": Total Sample (MinBias)");
	#endif 

	frames->Draw();

	//get the dataset with sWeights
	// The SPlot class adds a new variable that has the name of the corresponding yield + "_sw".
	RooDataSet* dataW = (RooDataSet*) ws.data("dataWithSWeights");
	RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
	RooDataSet* dataWBg = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_combinatorial_sw");

	RooPlot* frames2Bp = variable->frame();
	frames2Bp->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/nbin));
	frames2Bp->GetXaxis()->SetTitle(label);

	RooPlot* frames2Bg = variable->frame();
	frames2Bg->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/nbin));
	frames2Bg->GetXaxis()->SetTitle(label);
  
  	dataWBp->plotOn(frames2Bp, DataError(RooAbsData::SumW2),Binning(nbin));
  	dataWBg->plotOn(frames2Bg, DataError(RooAbsData::SumW2),Binning(nbin));

	#if Dataset == 0 // Data
	frames2Bp->SetTitle(label + " Signal-Data (splot)");
	frames2Bg->SetTitle(label + " Background-Data (splot)");
	#elif Dataset == 1 // MC
	frames2Bp->SetTitle(label + " Signal-MC (splot)");
	frames2Bg->SetTitle(label + " Background-MC (splot)");
	#elif Dataset == 2 // MinBias
	frames2Bp->SetTitle(label + " Signal-MinBias (splot)");
	frames2Bg->SetTitle(label + " Background-MinBias (splot)");
	#endif 
	 
	cdata->cd(3);  frames2Bp->Draw();
	cdata->cd(4);  frames2Bg->Draw();
	//------------------------------------
	// S A V E   P L O T S
	//------------------------------------
#if particle == 0 // Dstar
	#if Dataset == 0 // Data
	cdata->SaveAs("DstarData/SPlot"+label+"DATA.pdf");
	cdata->SaveAs("DstarData/SPlot"+label+"DATA.png");
	#elif Dataset == 1 // MC
	cdata->SaveAs("DstarMC/SPlot"+label+"MC.pdf");
	cdata->SaveAs("DstarMC/SPlot"+label+"MC.png");
	#elif Dataset == 2 // MinBias
	cdata->SaveAs("DstarMinBias/SPlot"+label+"MinBias.pdf");
	cdata->SaveAs("DstarMinBias/SPlot"+label+"MinBias.png");
	#endif
#elif particle == 1 // D0
	#if Dataset == 0 // Data
	cdata->SaveAs("D0DATA/SPlot"+label+"DATA.pdf");
	cdata->SaveAs("D0DATA/SPlot"+label+"DATA.png");
	#elif Dataset == 1 // MC
	cdata->SaveAs("D0MC/SPlot"+label+"MC.pdf");
	cdata->SaveAs("D0MC/SPlot"+label+"MC.png");
	#elif Dataset == 2 // MinBias
	cdata->SaveAs("D0MinBias/SPlot"+label+"MinBias.pdf");
	cdata->SaveAs("D0MinBias/SPlot"+label+"MinBias.png");
	#endif		
#endif
	//------------------------------------
	// H I S T O G R A M S
	//------------------------------------
#if particle == 0 // Dstar
	create_histogram(*Dsmass,"SPLOTDsmass", 1., dataWBg, dataWBp, data, 100, true);
	create_histogram(*variable,"SPLOT"+label, 1., dataWBg, dataWBp, data, 100, true);
#elif particle == 1 // D0
	create_histogram(*D0Kpimass,"SPLOTD0mass", 1., dataWBg, dataWBp, data, 100, true);
	create_histogram(*variable,"SPLOT"+label, 1., dataWBg, dataWBp, data, 100, true);	
#endif

	cout << "---------------------------------" << endl;
	cout << "End MakePlots" << endl;
	cout << "---------------------------------" << endl;
}
//=====================================================================
void build_pdf(RooWorkspace& w, std::string choice) {
	cout << "---------------------------------" << endl;
	cout << "Begin build_pdf()" << endl;
	cout << "---------------------------------" << endl;

#if particle == 0 // Dstar 
	RooRealVar Dmass = *(w.var("Dsmass"));
#elif particle == 1 // D0
	RooRealVar Dmass = *(w.var("D0Kpimass"));	
#endif
	RooDataSet* data = (RooDataSet*) w.data("data");
	double left = 1.99; 
	double right = 2.035;
	double mass_peak = 2.010 ;
	//----------------------------------
	//S I G N A L
	//----------------------------------
	RooRealVar mean("mean","mean",mass_peak,mass_peak-0.1,mass_peak+0.1);
	RooRealVar sigma1("sigma1","sigma1",0.0252,0.020,0.030);
	RooRealVar sigma2("sigma2","sigma2",0.01052,0.010,0.020);

	RooGaussian signal1("signal1","signal_gauss1",Dmass,mean,sigma1);
	RooGaussian signal2("signal2","signal_gauss2",Dmass,mean,sigma2);

	RooRealVar cofs("cofs", "cofs", 0.317, 0., 1.);
	RooAddPdf signal("signal", "signal", RooArgList(signal1,signal2),cofs);
	//sigma1.setConstant();
	//sigma2.setConstant();
	//cofs.setConstant();
	//---------------------------
	//systematics->signal
	//--------------------------
	RooRealVar sigma3("sigma3","sigma3",0.012,0.010,0.030);
	RooGaussian signal3("signal3","signal3",Dmass, mean, sigma3);
	//----------------------------
	//B A C K G R O U N D
	//---------------------------------------------
	//error function
	RooRealVar m_nonprompt_scale("m_nonprompt_scale", "m_nonprompt_scale", 4.74168e-02, 0, 1);
	RooRealVar m_nonprompt_shift("m_nonprompt_shift", "m_nonprompt_shift", 5.14425, 4.5, 6.); 
	m_nonprompt_shift.setConstant(kTRUE);
	m_nonprompt_scale.setConstant(kTRUE);

#if particle == 0 // Dstar
	RooGenericPdf erf("erf","erf", "TMath::Erfc((Dsmass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(Dmass, m_nonprompt_scale, m_nonprompt_shift));
#elif particle == 1 // D0
	RooGenericPdf erf("erf","erf", "TMath::Erfc((D0Kpimass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(Dmass, m_nonprompt_scale, m_nonprompt_shift));
#endif
	
	//exponential
	RooRealVar lambda("lambda","lambda",-2.,-5.,0.0);
	RooExponential fit_side("fit_side", "fit_side_exp", Dmass, lambda);

	//D meson component
	RooRealVar D_mean1("D_mean1","D_mean1",5.34693e+00,Dmass.getAsymErrorLo(),Dmass.getAsymErrorHi());
	RooRealVar D_mean2("D_mean2","D_mean2",5.46876e+00,Dmass.getAsymErrorLo(),Dmass.getAsymErrorHi());
	RooRealVar D_mean3("D_mean3","D_mean3",5.48073e+00,Dmass.getAsymErrorLo(),Dmass.getAsymErrorHi());

	RooRealVar sigma1_left("sigma1_left","sigma1_left",2.90762e-02,0.010,0.150);
	RooRealVar sigma1_right("sigma1_right","sigma1_right",6.52519e-02,0.010,0.150);

	RooRealVar Dsigma2("Dsigma2","Dsigma2",9.94712e-02,0.020,0.500);
	RooRealVar Dsigma3("Dsigma3","Dsigma3",3.30152e-01,0.020,0.500);

	RooRealVar fraction2("fraction2","fraction2",2.34646e-01,0.0,1.0);
	RooRealVar fraction3("fraction3","fraction3",1.14338e-01,0.0,1.0);

	D_mean1.setConstant(kTRUE);
	D_mean2.setConstant(kTRUE);
	D_mean3.setConstant(kTRUE);
	sigma1_left.setConstant(kTRUE);
	sigma1_right.setConstant(kTRUE);
	Dsigma2.setConstant(kTRUE);
	Dsigma3.setConstant(kTRUE);
	fraction2.setConstant(kTRUE);
	fraction3.setConstant(kTRUE);

	RooBifurGauss gaussian1("gaussian1","gaussian1",Dmass,D_mean1,sigma1_left,sigma1_right);
	RooGaussian gaussian2("gaussian2","gaussian2",Dmass,D_mean2,Dsigma2);
	RooGaussian gaussian3("gaussian3","gaussian3",Dmass,D_mean3,Dsigma3);
	RooAddPdf mesonD("mesonD","mesonD",RooArgList(gaussian3,gaussian2,
	gaussian1),RooArgList(fraction3,fraction2));

	Dmass.setRange("all", Dmass.getMin(),Dmass.getMax());
	Dmass.setRange("right",right,Dmass.getMax());
	Dmass.setRange("left",Dmass.getMin(),left);
	Dmass.setRange("peak",left,right);
	Dmass.setRange("peakright",left,Dmass.getMax());
	//n values
#if particle == 0 // Dstar
	double n_signal_initial = data->sumEntries(TString::Format("abs(Dsmass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(Dsmass-%g)<0.10&&abs(Dsmass-%g)>0.05",mass_peak,mass_peak));
#elif particle == 1 // D0
	double n_signal_initial = data->sumEntries(TString::Format("abs(D0Kpimass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(D0Kpimass-%g)<0.10&&abs(D0Kpimass-%g)>0.05",mass_peak,mass_peak));
#endif
	
	double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
	RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
	RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
	RooRealVar f_erf("f_erf","f_erf",2.50259e-01,0,1);
	RooProduct n_erf("n_erf","n_erf",RooArgList(n_signal,f_erf));

	RooRealVar f_mesonD("f_mesonD","f_mesonD",4.1E-5/1.026E-3,0.,0.1); 
	f_mesonD.setConstant(kTRUE);
	RooProduct n_mesonD("n_mesonD","n_mesonD",RooArgList(n_signal,f_mesonD));
	//systematics->bkg
	RooRealVar slope("slope","slope",0,-10,10);
	RooPolynomial poly_bkg("poly_bkg", "poly_bkg", Dmass, slope);

	//----------------------------------------------------------------
	//M O D E L   C H O I C E
	//----------------------------------------------------------------
	if(choice == "nominal"){
      RooAddPdf model("model", "model", RooArgList(signal,fit_side,erf,mesonD), RooArgList(n_signal,n_combinatorial,n_erf,n_mesonD));
		w.import(model); }
	else if(choice=="bkg_poly"){
		RooAddPdf model("model", "model", RooArgList(signal,poly_bkg,erf,mesonD), RooArgList(n_signal,n_combinatorial,n_erf,n_mesonD));
		w.import(model);}
	else if (choice=="bkg_range"){
		RooAddPdf model("model", "model", RooArgList(signal,fit_side, mesonD), RooArgList(n_signal,n_combinatorial, n_mesonD));
      w.import(model);}
	else if (choice == "signal1gauss"){
		RooAddPdf model("model", "model", RooArgList(signal3,fit_side,erf,mesonD), RooArgList(n_signal,n_combinatorial,n_erf,n_mesonD));
		w.import(model);}

	cout << "---------------------------------" << endl;
	cout << "End build_pdf()" << endl;
	cout << "---------------------------------" << endl;
}
//=====================================================================
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
#if particle == 0 // Dstar
	RooRealVar* Dspt = ws.var("Dspt");
  	RooRealVar Dsmass = *(ws.var("Dsmass"));
#elif particle == 1 // D0
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
#if particle == 0 // Dstar
		data_pt = (RooDataSet*) data->reduce(Form("Dspt>%lf",pt_bins[i]));
		data_pt = (RooDataSet*) data_pt->reduce(Form("Dspt<%lf",pt_bins[i+1]));
#elif particle == 1 // D0
		data_pt = (RooDataSet*) data->reduce(Form("D0Kpipt>%lf",pt_bins[i]));
		data_pt = (RooDataSet*) data_pt->reduce(Form("D0Kpipt<%lf",pt_bins[i+1]));
#endif
			
		ws.import(*data_pt, Rename(Form("data_pt_%d",i)));
   
		//perform fit and save result
		fit_pt = model->fitTo(*data_pt, Minos(true), Save());

		//plots the fit result
#if particle == 0 // Dstar
		RooPlot* massframe = Dsmass.frame(Title(""));
#elif particle == 1 // D0
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

#if particle == 0 // Dstar
	#if Dataset == 0 // Data
		b.SaveAs(Form("DstarData/PT%i_%iDstarDATA.pdf", int1, int2));
		b.SaveAs(Form("DstarData/PT%i_%iDstarDATA.png", int1, int2));
	#elif Dataset == 1 // MC
		b.SaveAs(Form("DstarMC/PT%i_%iDstarMC.pdf", int1, int2));
		b.SaveAs(Form("DstarMC/PT%i_%iDstarMC.png", int1, int2));
	#elif Dataset == 2 // MinBias
		b.SaveAs(Form("DstarMinBias/PT%i_%iDstarMinBias.pdf", int1, int2));
		b.SaveAs(Form("DstarMinBias/PT%i_%iDstarMinBias.png", int1, int2));
	#endif
#elif particle == 1 // D0
	#if Dataset == 0 // Data
		b.SaveAs(Form("D0DATA/PT%i_%iD0DATA.pdf", int1, int2));
		b.SaveAs(Form("D0DATA/PT%i_%iD0DATA.png", int1, int2));
	#elif Dataset == 1 // MC
		b.SaveAs(Form("D0MC/PT%i_%iD0MC.pdf", int1, int2));
		b.SaveAs(Form("D0MC/PT%i_%iD0MC.png", int1, int2));
	#elif Dataset == 2 // MinBias
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
#if particle == 0 // Dstar
			mean_w = data_wb->mean(*Dspt);
		   mean_s = data_pt->mean(*Dspt);
		   pt_mean[i] = data_wb->mean(*Dspt);
#elif particle == 1 // D0
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
	#if particle == 0 // Dstar
	RooRealVar* Dseta = ws.var("Dseta");
  	RooRealVar Dsmass = *(ws.var("Dsmass"));
	#elif particle == 1 // D0
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
		#if particle == 0 // Dstar
		data_eta = (RooDataSet*) data->reduce(Form("fabs(Dseta)>%lf",eta_bins[i]));
		data_eta = (RooDataSet*) data_eta->reduce(Form("fabs(Dseta)<%lf",eta_bins[i+1]));
		#elif particle == 1 // D0
		data_eta = (RooDataSet*) data->reduce(Form("fabs(D0Kpieta)>%lf",eta_bins[i]));
		data_eta = (RooDataSet*) data_eta->reduce(Form("fabs(D0Kpieta)<%lf",eta_bins[i+1]));
		#endif
		
		ws.import(*data_eta, Rename(Form("data_eta_%d",i)));
   
		//perform fit and save result
		fit_eta = model->fitTo(*data_eta, Minos(true), Save());

		//plots the fit result
		#if particle == 0 // Dstar
		RooPlot* massframe = Dsmass.frame(Title(""));
		#elif particle == 1 // D0
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
#if particle == 0 // Dstar
	#if Dataset == 0 // Data
		b.SaveAs(Form("DstarData/ETA%i_%iDstarDATA.png", lowbins, highbins));
		b.SaveAs(Form("DstarData/ETA%i_%iDstarDATA.pdf", lowbins, highbins));
	#elif Dataset == 1 // MC
		b.SaveAs(Form("DstarMC/ETA%i_%iDstarMC.png", lowbins, highbins));
		b.SaveAs(Form("DstarMC/ETA%i_%iDstarMC.pdf", lowbins, highbins));
	#elif Dataset == 2 // MinBias
		b.SaveAs(Form("DstarMinBias/ETA%i_%iDstarMinBias.png", lowbins, highbins));
		b.SaveAs(Form("DstarMinBias/ETA%i_%iDstarMinBias.pdf", lowbins, highbins));
	#endif
#elif particle == 1 // D0
	#if Dataset == 0 // Data
		b.SaveAs(Form("D0Data/ETA%i_%iD0DATA.png", lowbins, highbins));
		b.SaveAs(Form("D0Data/ETA%i_%iD0DATA.pdf", lowbins, highbins));
	#elif Dataset == 1 // MC
		b.SaveAs(Form("D0MC/ETA%i_%iD0MC.png", lowbins, highbins));
		b.SaveAs(Form("D0MC/ETA%i_%iD0MC.pdf", lowbins, highbins));
	#elif Dataset == 2 // MinBias
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

#if particle == 0 // Dstar
	RooRealVar Dslifetime = *(ws.var("Dslifetime"));
#elif particle == 1 // D0
	RooRealVar D0lifetime = *(ws.var("D0lifetime"));
#endif

	//RooDataSet* data_central = (RooDataSet*) ws.data("data_central");
	//Get the dataset with Wights
	RooDataSet* data = (RooDataSet*) ws.data("data");
	RooDataSet* dataW = (RooDataSet*) ws.data("dataWithSWeights");
	//Get the from the "dataWithSWeights" dataser only the signal
	RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
	RooDataSet* dataWBp2;

	RooDataSet* data_central = (RooDataSet*) ws.data("data_central");
	
#if particle == 0 // Dstar
	#if LifeSideBand == 0
	dataWBp2 = (RooDataSet*) dataWBp->reduce("Dslifetime>(0.0*pow(10,-12))");	
	#elif LifeSideBand == 1
	dataWBp2 = (RooDataSet*) data_central->reduce("Dslifetime>(0.0*pow(10,-12))");
	#endif
#elif particle == 1 // D0
	#if LifeSideBand == 0
	dataWBp2 = (RooDataSet*) dataWBp->reduce("D0lifetime>(0.0*pow(10,-12))");
	#elif LifeSideBand == 1
	dataWBp2 = (RooDataSet*) data_central->reduce("D0lifetime>(0.0*pow(10,-12))");
	#endif
#endif

	//------------------------------------------------
	// P D F    L I F E T I M E
	//------------------------------------------------
	Double_t lambdaPrompt = -1./(410.1*pow(10,-15)); //D0 decay   -1./(410.1*pow(10,-15))    -1./1.02525*pow(10,-13))
	RooRealVar lambda1("lambda1", "lambda1", lambdaPrompt);
	//RooRealVar lambda1( "lambda1", "lambda1", -1./(410.1*pow(10,-15)), -4./(410.1*pow(10,-15)), 4./(410.1*pow(10,-15)) );
	Double_t lambdaNonPrompt = -1./(1.5190*pow(10,-12)); //B decay
	RooRealVar lambda2("lambda2","lambda2", lambdaNonPrompt);	
#if particle == 0 // Dstar
	//RooExponential lifetime1("lifetime1", "lifetime1", Dslifetime, lambda1);
	RooExponential lifetime1("lifetime1", "lifetime1", Dslifetime, lambda1);
  	RooExponential lifetime2("lifetime2", "lifetime2", Dslifetime, lambda2);
#elif particle == 1 // D0
	RooExponential lifetime1("lifetime1", "lifetime1", D0lifetime, lambda1);
  	RooExponential lifetime2("lifetime2", "lifetime2", D0lifetime, lambda2);
#endif
  	
   lambda1.setConstant();
   lambda2.setConstant();
	RooRealVar prompt("prompt","prompt", 1000., 0.0, 5000000.);
	RooRealVar non_prompt("non_prompt","non_prompt", 1000., 0.0, 5000000.);
	RooAddPdf LifetimeModel("LifetimeModel","LifetimeModel", RooArgList(lifetime1,lifetime2), RooArgList(prompt,non_prompt));
	
	ws.import(LifetimeModel);

#if particle == 0 // Dstar
	#if Dataset == 0 // Data
	Dslifetime.setRange("range1",0.3*pow(10,-12),1.6*pow(10,-12));
	#elif Dataset == 1 // MC
	Dslifetime.setRange("range1",0.3*pow(10,-12),0.9*pow(10,-12));
	#elif Dataset == 2 // MinBias
	Dslifetime.setRange("range1",0.3*pow(10,-12),0.9*pow(10,-12));
	#endif
#elif particle == 1 // D0
	#if Dataset == 0 // Data
	D0lifetime.setRange("range1",0.8*pow(10,-12),5.0*pow(10,-12));
	#elif Dataset == 1 // MC
	D0lifetime.setRange("range1",0.6*pow(10,-12),2.8*pow(10,-12));
	#elif Dataset == 2 // MinBias
	D0lifetime.setRange("range1",0.7*pow(10,-12),5.0*pow(10,-12));
	#endif		
#endif
	//------------------------------------------------
	// F I T   D A T A
	//------------------------------------------------
	//LifetimeModel.fitTo(*dataWBp2,Range("all"));
	//LifetimeModel.fitTo(*dataWBp2,Extended());
	LifetimeModel.fitTo(*dataWBp2,Range("range1"));
	//------------------------------------------------
	// C A N V A S
	//------------------------------------------------
	TCanvas d;
  	d.SetTitle("");
#if particle == 0 // Dstar
	RooPlot* massframe = Dslifetime.frame();	
	#if LifeSideBand == 0
	massframe->SetTitle("D0(from D*) Lifetime - sPLOT");
	#elif LifeSideBand == 1
	massframe->SetTitle("D0(from D*) Lifetime - SideBand");
	#endif
#elif particle == 1 // D0
	RooPlot* massframe = D0lifetime.frame();
	#if LifeSideBand == 0	
	massframe->SetTitle("D0 Lifetime - sPLOT");
	#elif LifeSideBand == 1
	massframe->SetTitle("D0 Lifetime - SideBand");
	#endif
#endif

	//d.SetLogy();
	//massframe->GetYaxis()->SetLogy();
	//massframe->SetLogy();

	//------------------------------------------------
	// P L O T S 
	//------------------------------------------------
	dataWBp2->plotOn(massframe, RooFit::Name("dataWBp2") );
	double totalEntries = dataWBp2->sumEntries();
	LifetimeModel.plotOn(massframe, RooFit::Name("Fit"),Range("range1"),LineColor(kRed),LineStyle(1),LineWidth(2));
	double chis = massframe->chiSquare();
	LifetimeModel.plotOn(massframe, RooFit::Name("lifetime1"),Components("lifetime1"),Range("range1"),LineColor(kBlue),LineStyle(kDashed));
	LifetimeModel.plotOn(massframe, RooFit::Name("lifetime2"),Components("lifetime2"),Range("range1"),LineColor(kOrange),LineStyle(kDashed));
	LifetimeModel.paramOn(massframe,Layout(0.75,0.99,0.75));
	massframe->getAttText()->SetTextSize(0.028);
	massframe->GetYaxis()->SetTitleOffset(1.3);
	massframe->SetXTitle("Lifetime [s]");
	massframe->Draw();

	//------------------------------------------------
	// L E G E N D S
	//------------------------------------------------
	TLatex* tex13 = new TLatex(0.15, 0.85, Form("#chi^{2}/DOF = %.3lf",chis));
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04);	
	tex13->Draw();
	//TLegend *leg = new TLegend (0.4, 0.5, 0.6, 0.7); 
	TLegend *leg = new TLegend (0.45,0.70,0.65,0.85);
	#if Dataset == 0 // Data
	leg->AddEntry(massframe->findObject("dataWBp2"), "Data", "LP");
	#elif Dataset == 1 // MC
	leg->AddEntry(massframe->findObject("dataWBp2"), "MC", "LP");
	#elif Dataset == 2 // MinBias
	leg->AddEntry(massframe->findObject("dataWBp2"), "MinBias", "LP");
	#endif
	leg->AddEntry("dataWBp2",(Form("Entries: %2.0f", totalEntries))," ");
	leg->AddEntry("Fit","Fit","L");
	leg->AddEntry("lifetime1","Prompt","L");
	leg->AddEntry("lifetime2","Non-Prompt","L");
	leg->Draw();
	//------------------------------------------------
	// S A V E   C A N V A S
	//------------------------------------------------
#if particle == 0 // Dstar
	#if Dataset == 0 // Data
	d.SaveAs("DstarData/DsLifetimeDATA.pdf");
	d.SaveAs("DstarData/DsLifetimeDATA.png");
	#elif Dataset == 1 // MC
	d.SaveAs("DstarMC/DsLifetimeMC.pdf");
	d.SaveAs("DstarMC/DsLifetimeMC.png");
	#elif Dataset == 2 // MinBias
	d.SaveAs("DstarMinBias/DsLifetimeMinBias.pdf");
	d.SaveAs("DstarMinBias/DsLifetimeMinBias.png");
	#endif
#elif particle == 1 // D0
	#if Dataset == 0 // Data
	d.SaveAs("D0DATA/D0LifetimeDATA.pdf");
	d.SaveAs("D0DATA/D0LifetimeDATA.png");
	#elif Dataset == 1 // MC
	d.SaveAs("D0MC/D0LifetimeMC.pdf");
	d.SaveAs("D0MC/D0LifetimeMC.png");
	#elif Dataset == 2 // MinBias
	d.SaveAs("D0MinBias/D0LifetimeMinBias.pdf");
	d.SaveAs("D0MinBias/D0LifetimeMinBias.png");
	#endif
#endif
	
	cout << "lambda1: "<< lambda1.getVal() << endl;
	cout << "tau1: "<< (-1./lambda1.getVal()) << " s" << endl;
	cout << "---------------------------------" << endl;
	cout << "End ForLifetime" << endl;
	cout << "---------------------------------" << endl;	
}
//=========================================================================
void GenLifetime(RooWorkspace& ws){
	cout << "---------------------------------" << endl;
	cout << "GenLifetime" << endl;
	cout << "---------------------------------" << endl;
	
#if particle == 0 // Dstar
	//RooRealVar MCD0lifetime = *(ws.var("MCD0lifetime"));
	RooRealVar MCD0lifetimeMatching = *(ws.var("MCD0lifetimeMatching"));
	RooDataSet* data2 = (RooDataSet*) ws.data("data2");
	//data2 = (RooDataSet*) data2->reduce("MCD0lifetime>(0.0*pow(10,-12))"); 
	data2 = (RooDataSet*) data2->reduce("MCD0lifetimeMatching>(0.0*pow(10,-12))");
#elif particle == 1 // D0
	//RooRealVar MCpromptD0lifetime = *(ws.var("MCpromptD0lifetime"));
	RooRealVar MCpromptD0lifetimeMatching = *(ws.var("MCpromptD0lifetimeMatching"));
	RooDataSet* data2 = (RooDataSet*) ws.data("data2");
	//data2 = (RooDataSet*) data2->reduce("MCpromptD0lifetime>(0.0*pow(10,-12))");
	data2 = (RooDataSet*) data2->reduce("MCpromptD0lifetimeMatching>(0.0*pow(10,-12))");
#endif
	
	//------------------------------------------------
	// P D F    L I F E T I M E
	//------------------------------------------------
	Double_t lambdaPrompt = -1./(410.1*pow(10,-15)); //D0 decay
	RooRealVar lambda1("lambda1", "lambda1", lambdaPrompt);
	Double_t lambdaNonPrompt = -1./(1.5190*pow(10,-12)); //B decay
	RooRealVar lambda2("lambda2","lambda2", lambdaNonPrompt);
#if particle == 0 // Dstar
	//RooExponential lifetime1("lifetime1", "lifetime1", MCD0lifetime, lambda1);
	RooExponential lifetime1("lifetime1", "lifetime1", MCD0lifetimeMatching, lambda1);
#elif particle == 1 // D0
	//RooExponential lifetime1("lifetime1", "lifetime1", MCpromptD0lifetime, lambda1);
	RooExponential lifetime1("lifetime1", "lifetime1", MCpromptD0lifetimeMatching, lambda1);
#endif
	
   lambda1.setConstant();

	RooRealVar prompt("prompt","prompt", 1000., 0.0, 5000000.);
	RooRealVar non_prompt("non_prompt","non_prompt", 1000., 0.0, 5000000.);
	RooAddPdf LifetimeModel("LifetimeModel","LifetimeModel", RooArgList(lifetime1), RooArgList(prompt));
	ws.import(LifetimeModel);
#if particle == 0 // Dstar
	//MCD0lifetime.setRange("range1",0.0*pow(10,-12),1.6*pow(10,-12));
	MCD0lifetimeMatching.setRange("range1",0.3*pow(10,-12),1.6*pow(10,-12));
#elif particle == 1 // D0
	//MCpromptD0lifetime.setRange("range1",0.0*pow(10,-12),1.6*pow(10,-12));
	MCpromptD0lifetimeMatching.setRange("range1",0.6*pow(10,-12),2.6*pow(10,-12));
#endif

	//------------------------------------------------
	// F I T   D A T A
	//------------------------------------------------
	LifetimeModel.fitTo(*data2,Range("range1"));	
	//------------------------------------------------
	// C A N V A S
	//------------------------------------------------
	TCanvas d;
  	d.SetTitle("");
#if particle == 0 // Dstar
	//RooPlot* massframe = MCD0lifetime.frame();
	RooPlot* massframe = MCD0lifetimeMatching.frame();	
	massframe->SetTitle("Matched Gen D0(from D*) Lifetime");
#elif particle == 1 // D0
	//RooPlot* massframe = MCpromptD0lifetime.frame();
	RooPlot* massframe = MCpromptD0lifetimeMatching.frame();
	massframe->SetTitle("Matched Gen D0 Lifetime");
#endif
	
	//d.SetLogy();
	//massframe->GetYaxis()->SetLogy();
	//massframe->SetLogy();
	//------------------------------------------------
	// P L O T S 
	//------------------------------------------------
	data2->plotOn(massframe, RooFit::Name("data2") );
	double totalEntries = data2->sumEntries();
	LifetimeModel.plotOn(massframe, RooFit::Name("prompt"),Range("range1"),LineColor(kRed),LineStyle(1),LineWidth(2));
	double chis = massframe->chiSquare();
	LifetimeModel.paramOn(massframe,Layout(0.75,0.99,0.75));
	massframe->getAttText()->SetTextSize(0.028);
	massframe->GetYaxis()->SetTitleOffset(1.3);
	massframe->SetXTitle("Lifetime [s]");
	massframe->Draw();
	//------------------------------------------------
	// L E G E N D S
	//------------------------------------------------
	TLatex* tex13 = new TLatex(0.15, 0.85, Form("#chi^{2}/DOF = %.3lf",chis));
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04);	
	tex13->Draw();
	//TLegend *leg = new TLegend (0.4, 0.5, 0.6, 0.7); 
	TLegend *leg = new TLegend (0.45,0.70,0.65,0.85);
	leg->AddEntry(massframe->findObject("data2"), "Gen", "LP");
	leg->AddEntry("data2",(Form("Entries: %2.0f", totalEntries))," ");
	leg->AddEntry("prompt","prompt","L");
	leg->Draw();
	//------------------------------------------------
	// S A V E   C A N V A S
	//------------------------------------------------
#if particle == 0 // Dstar
	d.SaveAs("DstarMC/DsLifeTimeGen.pdf");
	d.SaveAs("DstarMC/DsLifeTimeGen.png");
#elif particle == 1 // D0
	d.SaveAs("D0MC/D0LifeTimeGen.pdf");
	d.SaveAs("D0MC/D0LifeTimeGen.png");
#endif
	cout << "---------------------------------" << endl;
	cout << "End GenLifetime" << endl;
	cout << "---------------------------------" << endl;
}
//=========================================================================
void ForLifetimeV2(RooWorkspace& ws){
//Correctly way to fit Proper Decay after the selection criteria
	cout << "---------------------------------" << endl;
	cout << "ForLifetimeV2" << endl;
	cout << "---------------------------------" << endl;

	RooRealVar DstarCt = *(ws.var("DstarCt"));
	//RooDataSet* data_central = (RooDataSet*) ws.data("data_central");
	//Get the dataset with Wights
	RooDataSet* data = (RooDataSet*) ws.data("data");
	RooDataSet* dataW = (RooDataSet*) ws.data("dataWithSWeights");
	//Get the from the "dataWithSWeights" dataser only the signal
	RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");

	//------------------------------------------------
	// P D F 
	//------------------------------------------------
	//RooRealVar DstarCt("DstarCt", "DstarCt", 0, 1.6);
   RooRealVar cterr("cterr", "per-event error on ct", 0.04);
   // Build a gaussian resolution model scaled by the per-event error = gauss(ct,bias,sigma*cterr)
   RooRealVar bias("bias", "bias", 5, -10, 20);
   RooRealVar sigma("sigma", "per-event error scale factor", 0.5, 0.01, 1);
   RooGaussModel gm("gm1", "gauss model scaled bt per-event error", DstarCt, bias, sigma, cterr);

   // Construct model as  exponential(ct) (x) gauss(ct|cterr)
   RooRealVar tau("tau", "tau", 0.41, 0.001, 1.);
   RooDecay LifetimeModel("decay_gm", "decay", DstarCt, tau, gm, RooDecay::SingleSided);

	//------------------------------------------------
	// F I T   D A T A
	//------------------------------------------------
	LifetimeModel.fitTo(*dataWBp,Range("all"));
	//------------------------------------------------
	// C A N V A S
	//------------------------------------------------
	TCanvas d;
  	d.SetTitle("");

	RooPlot* massframe = DstarCt.frame();	
	massframe->SetTitle("D0(from D*) proper decay length - sPLOT");
	
	//------------------------------------------------
	// P L O T S 
	//------------------------------------------------
	dataWBp->plotOn(massframe, RooFit::Name("dataWBp") );
	double totalEntries = dataWBp->sumEntries();
	LifetimeModel.plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
	double chis = massframe->chiSquare();
	//LifetimeModel.plotOn(massframe, RooFit::Name("lifetime1"),Components("lifetime1"),Range("range1"),LineColor(kBlue),LineStyle(kDashed));
	//LifetimeModel.plotOn(massframe, RooFit::Name("lifetime2"),Components("lifetime2"),Range("range1"),LineColor(kOrange),LineStyle(kDashed));
	LifetimeModel.paramOn(massframe,Layout(0.75,0.99,0.75));
	massframe->getAttText()->SetTextSize(0.028);
	massframe->GetYaxis()->SetTitleOffset(1.3);
	massframe->SetXTitle("ct [cm]");
	massframe->Draw();

	//------------------------------------------------
	// L E G E N D S
	//------------------------------------------------
	TLatex* tex13 = new TLatex(0.15, 0.85, Form("#chi^{2}/DOF = %.3lf",chis));
	tex13->SetNDC(kTRUE); tex13->SetTextFont(42); tex13->SetTextSize(0.04);	
	tex13->Draw();
	//TLegend *leg = new TLegend (0.4, 0.5, 0.6, 0.7); 
	TLegend *leg = new TLegend (0.45,0.70,0.65,0.85);
	#if Dataset == 0 // Data
	leg->AddEntry(massframe->findObject("dataWBp"), "Data", "LP");
	#elif Dataset == 1 // MC
	leg->AddEntry(massframe->findObject("dataWBp"), "MC", "LP");
	#elif Dataset == 2 // MinBias
	leg->AddEntry(massframe->findObject("dataWBp"), "MinBias", "LP");
	#endif
	leg->AddEntry("dataWBp",(Form("Entries: %2.0f", totalEntries))," ");
	leg->AddEntry("Fit","Fit","L");
	//leg->AddEntry("lifetime1","Prompt","L");
	//leg->AddEntry("lifetime2","Non-Prompt","L");
	leg->Draw();

	//------------------------------------------------
	// S A V E   C A N V A S
	//------------------------------------------------
	#if Dataset == 0 // Data
	d.SaveAs("DstarData/DstarProperTimeDATA.pdf");
	d.SaveAs("DstarData/DstarProperTimeDATA.png");
	#elif Dataset == 1 // MC
	d.SaveAs("DstarMC/DstarProperTimeMC.pdf");
	d.SaveAs("DstarMC/DstarProperTimeMC.png");
	#elif Dataset == 2 // MinBias
	d.SaveAs("DstarMinBias/DstarProperTimeMinBias.pdf");
	d.SaveAs("DstarMinBias/DstarProperTimeMinBias.png");
	#endif
	
	cout << "---------------------------------" << endl;
	cout << "End ForLifetimeV2" << endl;
	cout << "---------------------------------" << endl;	
}
