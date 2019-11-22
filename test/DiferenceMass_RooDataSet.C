//PDF - Probability density function
#include "RooDataSet.h"
using namespace RooFit ;

void DiferenceMass_RooDataSet()
{
//OPENING A ROOT FILE AND GETTING A TH1D WICH CONTAINS THE DATA
TFile* file = new TFile("DMesonsHistograms.root","read"); 
//TH1D *massHist = (TH1D *)file->Get("DsMinusD0Histo");
cout << "Reading root file" << endl;

TTree *tree = (TTree *)file->Get("t_analysis");
//TBranch *branch = (TBranch *)tree->Get("D0mass");

RooRealVar DsMinusD0("DsMinusD0","DsMinusD0",0.14,0.16);
//RooDataHist rooMassHisto("rooMassHisto","rooMassHisto",mass,Import(*massHist));
RooDataSet ds("ds", "ds", DsMinusD0, Import(*tree));
cout << "reading histogram" << endl;


// -----------------------------------------
// I m p o r t   T T r e e   i n t o   a   R o o D a t a S e t
// -----------------------------------------------------------
// Define observable y
//RooRealVar y("y", "y", -10, 10);
//RooDataSet ds("ds", "ds", RooArgSet(x, y), Import(*tree));

//--------------------------------------------------
//R O O F I T   V A R I A B L E S
//--------------------------------------------------
RooRealVar mean1( "mean1", "mean1", 0.1354, 0.125, 0.1455);
RooRealVar sigma1( "sigma1", "sigma1", 0.00022, 0.00004, 0.00055);
RooRealVar sigma2( "sigma2", "sigma2", 0.0008, 0.00009, 0.0009);
cout << "ROOFIT VARIABLES OK" << endl;

//--------------------------------------------------
//F U N C T I O N S
//--------------------------------------------------
//Gaussian (S I G N A L)
RooGaussian gauss1("gauss1","gauss1",DsMinusD0,mean1,sigma1);
RooGaussian gauss2("gauss2","gauss2",DsMinusD0,mean1,sigma2);

//------------------------------------------------
// A D D  S I G N A L   A N D   B A C K G R O U N D
// ------------------------------------------------
// Sum the composite signal and background 
RooRealVar gauss1_fac("gauss1_fac","fraction of Gaussian1", 4000, 1100, 13500);
RooRealVar gauss2_fac("gauss2_fac","fraction of Gaussian2", 800, 600, 2000);
RooAddPdf  sig("sig","gauss1 + gauss2",RooArgList(gauss1,gauss2),RooArgList(gauss1_fac,gauss2_fac));
cout << "ROOFIT FUNCTIONS OK" << endl;

//------------------------------------------------
//F I T T I N G
//------------------------------------------------
// Print unbinned dataset with default frame binning (100 bins)
RooPlot *frame = DsMinusD0.frame(Title("#Deltam = m(K#pi#pi_{slow}) - m(K#pi)"));
frame->SetXTitle("#Delta m [GeV/c2]");
sig.fitTo(ds);
    
ds.plotOn(frame,Name("Data"));
sig.plotOn(frame,Name("SPLUSBK"),LineColor(kBlue));


RooArgSet * pars = sig.getParameters(ds);

int nfloatpars = pars->selectByAttrib("Constant",kFALSE)->getSize(); 
double mychsq = frame->chiSquare("SPLUSBK","Data", nfloatpars);

sig.paramOn(frame,Layout(0.55,0.99,0.85),Format("NE"),Label(Form("#chi^{2}/ndf = %f", mychsq)) );
frame->getAttText()->SetTextSize(0.03);//change the font size


cout << "--------------------------------------" << endl;
cout << "nfloatpars: "<< nfloatpars << endl;
cout << "mychsq: " << mychsq << endl;
cout << "--------------------------------------" << endl;


// Draw all frames on a canvas
TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 800);
frame->Draw();

double totalentries = ds.sumEntries();
cout << "totalentries: "<< totalentries << endl;

TLegend *leg1 = new TLegend(0.50,0.25,0.70,0.35);
leg1->SetFillColor(kWhite);
leg1->SetLineColor(kBlack);
leg1->AddEntry("Data","MC","LP");
leg1->AddEntry("Data",(Form("Entries: %2.0f", totalentries)),"LP");
leg1->Draw();

//----------------------------------------
//S A V E   C A N V A S   A S   P N G
//----------------------------------------
ds.Print();
//canvas->SaveAs("/eos/user/r/ragomesd/analysisB2019/rootfit/FitDsMinusD0_MC.png");
canvas->SaveAs("/home/raphael/cernbox/analysisB2019/rootfit/FitDsMinusD0_MC.png");
    


cout << "---------------------------------" << endl;
cout << "E N D   P R O G R A M" << endl;
cout << "---------------------------------" << endl;
}//end program
