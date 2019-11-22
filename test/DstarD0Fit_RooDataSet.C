//PDF - Probability density function
using namespace RooFit ;

void DstarD0Fit_RooDataSet()
{

//OPENING A ROOT FILE AND GETTING A TH1D WICH CONTAINS THE DATA
TFile* file = new TFile("DMesonsHistograms.root","read"); 
//TH1D *massHist = (TH1D *)file->Get("DsMinusD0Histo");
cout << "Reading root file" << endl;

TTree *tree = (TTree *)file->Get("t_analysis");
//TBranch *branch = (TBranch *)tree->Get("D0mass");

RooRealVar D0mass("D0mass","D0mass",1.77,1.95);
//RooRealVar mass("mass","mass",1.93,2.10);
RooRealVar D0pt("D0pt","D0pt",3.,20);
RooDataSet ds("ds", "ds", RooArgSet(D0mass, D0pt), Import(*tree));


cout << "reading histogram" << endl;

//--------------------------------------------------
//Roofit Variables
//--------------------------------------------------
RooRealVar mean1( "mean1", "mean1", 1.864, 1.840, 1.890);
//RooRealVar mean1( "mean1", "mean1", 2.010, 2.000, 2.100);
RooRealVar sigma1( "sigma1", "sigma1", 0.020, 0.016, 0.06);
RooRealVar sigma2( "sigma2", "sigma2", 0.010, 0.006, 0.04);
RooRealVar poly_c0("poly_c0", "coefficient of x^0 term", 290253, 20000, 340000);
RooRealVar poly_c1("poly_c1", "coefficient of x^1 term", -208009, -218000, -188000);
RooRealVar poly_c2("poly_c2", "coefficient of x^2 term", 50000, 30000, 80000);
cout << "Roofit variables ok" << endl;

//---------------------------------------------------------------
//F U N C T I O N S
//---------------------------------------------------------------
//GAUSSIAN (S I G N A L)
RooGaussian gauss1("gauss1","gauss1",D0mass,mean1,sigma1);
RooGaussian gauss2("gauss2","gauss2",D0mass,mean1,sigma2);
//ROOFIT third degree polynomial (B A C K G R O U N D)
RooPolynomial poly("poly", "poly", D0mass, RooArgList(poly_c0, poly_c1, poly_c2));

//------------------------------------------------
// A d d  s i g n a l   a n d   b a c k g r o u n d
// ------------------------------------------------
// SUM THE COMPOSITE SIGNAL AND BACKGROUND 
RooRealVar gauss1_frac("gauss1_frac","fraction of Gaussian1", 7000, 800, 8500);
RooRealVar gauss2_frac("gauss2_frac","fraction of Gaussian2", 7000, 800, 8500);
// BACKGROUND 
RooRealVar poly_frac("poly_frac","fraction of poly", 12300, 500 , 12600);
//RooAddPdf  gauss1signal("gauss1signal","gauss1signal",RooArgList(gauss1));
//RooAddPdf  sigbkg("sigbkg","gauss1 + poly",RooArgList(gauss1,poly),RooArgList(gauss1_frac,poly_frac));
//RooAddPdf  sig("sig","gauss1",RooArgList(gauss1),RooArgList(gauss1_frac));
RooAddPdf  sig("sig","gauss1 + gauss2",RooArgList(gauss1,gauss2),RooArgList(gauss1_frac,gauss2_frac));
cout << "Roofit Function ok" << endl;

//------------------------------------------------
//F I T T I N G
//------------------------------------------------
// Print unbinned dataset with default frame binning (100 bins)
RooPlot *frame = D0mass.frame(Title("Fit of the D0"));
frame->SetXTitle("Invariant Mass [GeV/c2]");
sig.fitTo(ds);

// Print unbinned dataset with default frame binning (100 bins)
ds.plotOn(frame,Name("Data")); ds.plotOn(frame, Binning(100));
sig.plotOn(frame,Name("SPLUSBK"),LineColor(kBlue));


RooArgSet * pars = sig.getParameters(ds);

int nfloatpars = pars->selectByAttrib("Constant",kFALSE)->getSize(); 
double mychsq = frame->chiSquare("SPLUSBK","Data", nfloatpars);

sig.paramOn(frame,Layout(0.65,0.99,0.85),Format("NE"),Label(Form("#chi^{2}/ndf = %f", mychsq)) );
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

TLegend *leg1 = new TLegend(0.15,0.75,0.35,0.85);
leg1->SetFillColor(kWhite);
leg1->SetLineColor(kBlack);
leg1->AddEntry("Data","MC","LP");
leg1->AddEntry("Data",(Form("Entries: %2.0f", totalentries)),"LP");
leg1->Draw();

//----------------------------------------
//S A V E   C A N V A S   A S   P N G
//----------------------------------------
ds.Print();
//canvas->SaveAs("/eos/user/r/ragomesd/analysisB2019/rootfit/FitD0fromDstar_MC.png");
canvas->SaveAs("/home/raphael/cernbox/analysisB2019/rootfit/FitD0fromDstar_MC.png");


cout << "---------------------------------" << endl;
cout << "E N D   P R O G R A M" << endl;
cout << "---------------------------------" << endl;
}//end program
