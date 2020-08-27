#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit ;

TH1* makeTH1() ;
TTree* makeTTree() ;

void test()
{
    TFile *file0  = TFile::Open("Run2011_upsilon.root");
    TTree *DataTree = (TTree*)file0->Get(("tagandprobe/AnalysisTree"));
    
    double _mmin = 9;  double _mmax = 11;

    RooRealVar InvariantMass("InvariantMass", "InvariantMass", _mmin, _mmax);
    RooDataSet *Data = new RooDataSet("DATA", "DATA", DataTree, RooArgSet(InvariantMass));
    
    TH1* hist_all = Data->createHistogram("Invariant_mass", InvariantMass);
    hist_all->SetAxisRange(_mmin, _mmax);
    hist_all->Draw();
    
    double mass_peak1 = 9.46030;
    double mass_peak2 = 10.02326;
    double mass_peak3 = 10.3552;
    
    TCanvas* c_all  = new TCanvas;
    
    // Create a binned dataset that imports contents of TH1 and associates itscontents to observable 'mass'
    RooDataHist dh("dh","dh",InvariantMass,RooFit::Import(*hist_all));
    //RooDataHist dh_pass("dh_pass","dh_pass",InvariantMass,RooFit::Import(*hist_pass));
    
    // O BACKGROUND AGORA ESTÁ MAL FITADO
    RooRealVar a0("a0", "a0", -1000, 1000);
    RooRealVar a1("a1", "a1", -1000, 1000);
    RooRealVar a2("a2", "a2", -1000, 1000);
    RooRealVar a3("a0", "a0", -1000, 1000);
    RooRealVar a4("a2", "a2", -1000, 1000);
    RooRealVar a5("a0", "a0", -1000, 1000);
    RooRealVar a6("a2", "a2", -1000, 1000);
    RooRealVar a7("a0", "a0", -1000, 1000);
    
    RooChebychev background("cpol","cpol", InvariantMass,RooArgList(a0,a1,a2, a3, a5, a6, a7));
    
    RooRealVar sigma("sigma","sigma",0.05*(_mmax-_mmin),0.,0.5*(_mmax-_mmin));
    
    RooRealVar mean1("mean1","mean1",_mmin,(mass_peak1+mass_peak2)/2.);
    RooRealVar mean2("mean2","mean2",(mass_peak1+mass_peak2)/2.,(mass_peak3+mass_peak2)/2.);
    RooRealVar mean3("mean3","mean3",(mass_peak3+mass_peak2)/2.,_mmax);
    //FIT FUNCTIONS
    
    // --Gaussian as the signal pdf
    RooGaussian gaussian1("signal1","signal1",InvariantMass,mean1,sigma);
    RooGaussian gaussian2("signal2","signal2",InvariantMass,mean2,sigma);
    RooGaussian gaussian3("signal3","signal3",InvariantMass,mean3,sigma);
    
    double n_signal_initial1 =(dh.sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak1)) -dh.sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak1,mass_peak1))) / dh.sumEntries();
    double n_signal_initial2 =(dh.sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak2)) -dh.sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak2,mass_peak2))) / dh.sumEntries();
    double n_signal_initial3 =(dh.sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak3)) -dh.sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak3,mass_peak3))) / dh.sumEntries();
    
    double n_signal_initial_total = n_signal_initial1 + n_signal_initial2 + n_signal_initial3;
    
    //signal = gaussian1*frac1 + gaussian2*frac2 + gaussian3*(1-(frac1 + frac2))
    //S(signal)d mass = 1
    RooRealVar frac1("frac1","frac1",0.333,-1.,1.);
    RooRealVar frac2("frac2","frac2",0.333,-1.,1.);
    //para RooArgList N-1, assume como frações
    RooAddPdf* signal;
    signal = new RooAddPdf("signal", "signal", RooArgList(gaussian1, gaussian2,gaussian3), RooArgList(frac1, frac2));
    
    double n_back_initial = 1. - n_signal_initial1 - n_signal_initial2 -n_signal_initial3;
    
    RooRealVar n_signal_total("n_signal_total","n_signal_total",n_signal_initial_total,0.,dh.sumEntries());
    
    RooRealVar n_back("n_back","n_back",n_back_initial,0.,dh.sumEntries());
    //para RooArgList N, assume como normalizações
    //modelo_total = n_signal_total*signal + n_back*background
    RooAddPdf* model;
    model = new RooAddPdf("model","model", RooArgList(*signal, background),RooArgList(n_signal_total, n_back));
    
    RooFitResult* fitres = new RooFitResult; //saves fit result
    fitres = model->fitTo(dh, RooFit::Save());
    
    RooRealVar* pass_mean1 = (RooRealVar*) fitres->floatParsFinal().find("mean1");
    RooRealVar* pass_mean2 = (RooRealVar*) fitres->floatParsFinal().find("mean2");
    RooRealVar* pass_mean3 = (RooRealVar*) fitres->floatParsFinal().find("mean3");
    
    RooRealVar* pass_sigma = (RooRealVar*) fitres->floatParsFinal().find("sigma");
    
    double mean1_value = pass_mean1->getVal();
    double mean2_value = pass_mean2->getVal();
    double mean3_value = pass_mean3->getVal();
    
    double sigma_value = pass_sigma->getVal();
    
    RooPlot *frame = InvariantMass.frame(RooFit::Title("Invariant Mass"));
    
    c_all->Draw();
    
    frame->SetTitle("ALL");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    frame->SetYTitle(Form("Events / %3.1f MeV/c^{2}",hist_all->GetBinWidth(1)*1000));
    dh.plotOn(frame);
            
    model->plotOn(frame);
    model->plotOn(frame,RooFit::Components("signal1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame,RooFit::Components("signal2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame,RooFit::Components("signal3"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
    model->plotOn(frame,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    //model->paramOn(frame,RooFit::Layout(0.55, 0.99, 0.8));
    frame->Draw("");
    
    c_all->SaveAs("fit.pdf");
}
