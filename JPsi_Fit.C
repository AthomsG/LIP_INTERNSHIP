#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCbShape.h"
#include "TEfficiency.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

using namespace RooFit;

//EFFICIENCY TEST
char* strcat(string destination, string source)
{
    int size = destination.size() + source.size();
    char* output = new char[size + 1];
    for (int i = 0; i < size; i++)
    {
        if (i < destination.size())
            output[i] = destination[i];
        else
            output[i] = source[i - destination.size()];
    }
    output[size] = '\0';
    return output;
}

double* doFit(string condition)
{
    TFile* file0            = new TFile("DATA/JPsi/TP_JPSI_DATA.root"); //PATH TO HISTOGRAM
    TTree *DataTree = (TTree*)file0->Get(("JPSI_DATA"));
    
    string file     = "Result/";
    string all_pdf  = "_ALL.pdf";
    string pass_pdf = "_PASS.pdf";
       
    RooRealVar PassingProbeStandAloneMuon("PassingProbeStandAloneMuon", "PassingProbeStandAloneMuon", 0, 1); //Muon_Id
       
    double _mmin = 2.8;  double _mmax = 3.3;
       
    RooRealVar InvariantMass("InvariantMass", "InvariantMass", _mmin, _mmax);
    RooRealVar ProbeMuon_Pt("ProbeMuon_Pt", "ProbeMuon_Pt", 0, 60);
    RooRealVar ProbeMuon_Eta("ProbeMuon_Eta", "ProbeMuon_Eta", -3, 3);
    RooRealVar ProbeMuon_Phi("ProbeMuon_Phi", "ProbeMuon_Phi", -2, 2);
       
    RooFormulaVar* redeuce = new RooFormulaVar("PPTM", condition.c_str(), RooArgList(ProbeMuon_Pt));
    RooDataSet *Data_ALL    = new RooDataSet("DATA", "DATA", DataTree, RooArgSet(InvariantMass, PassingProbeStandAloneMuon, ProbeMuon_Pt), *redeuce);
    RooFormulaVar* cutvar = new RooFormulaVar("PPTM", strcat(condition.c_str(), "&& PassingProbeStandAloneMuon == 1"), RooArgList(PassingProbeStandAloneMuon, ProbeMuon_Pt));
    RooDataSet *Data_PASSING = new RooDataSet("DATA_PASS", "DATA_PASS", DataTree, RooArgSet(InvariantMass, PassingProbeStandAloneMuon, ProbeMuon_Pt), *cutvar);//
    
    //FITTING PARAMETERS
       
    double mass_peak1 = 3.1;
       
    TCanvas* c_all  = new TCanvas;
    TCanvas* c_pass = new TCanvas;
       
    //DECLARE OBSERVABLE X - INVARIANT MASS BETWEEN _mmin AND _mmax
       
    // Create a binned dataset that imports contents of TH1 and associates itscontents to observable 'mass'
    RooDataHist* dh      = Data_ALL->binnedClone();
    RooDataHist* dh_pass = Data_PASSING->binnedClone();

    RooRealVar lambda("lambda","lambda",-1.91, -2.,-1);
    RooExponential background("background", "background", InvariantMass, lambda);
       
    RooRealVar mean("mean","mean",3.094,3.085,3.1);
    RooRealVar sigma_cb("sigma_cb","sigma_cb", 0.038, 0.035, 0.04);
    RooRealVar alpha("alpha", "alpha", 1.71, 1.7, 1.8);
    RooRealVar n("n", "n", 3.96, 3.9, 4.);
    n.setConstant(kTRUE);
       
    //FIT FUNCTIONS
    RooRealVar sigma("sigma","sigma",0.05*(_mmax-_mmin),0.,0.5*(_mmax-_mmin));
    RooGaussian gaussian("GS","GS",InvariantMass,mean,sigma);
    RooCBShape crystalball("CB", "CB", InvariantMass, mean, sigma_cb, alpha, n);
       
    double n_signal_initial_total =10000;
       
    //signal = gaussian1*frac1 + gaussian2*frac2 + gaussian3*(1-(frac1 + frac2))
    //S(signal)d mass = 1
    RooRealVar frac("frac","frac",0.5,0.,1.);
    //para RooArgList N-1, assume como frações
    RooAddPdf* signal;
    signal = new RooAddPdf("signal", "signal", RooArgList(gaussian, crystalball), RooArgList(frac));
       
    double n_back_initial = 1. - n_signal_initial_total;
       
    RooRealVar n_signal_total("n_signal_total","n_signal_total",n_signal_initial_total,0.,dh->sumEntries());
       
    RooRealVar n_back("n_back","n_back",n_back_initial,0.,dh->sumEntries());
    //para RooArgList N, assume como normalizações
    //modelo_total = n_signal_total*signal + n_back*background
    RooAddPdf* model;
    model = new RooAddPdf("model","model", RooArgList(*signal, background),RooArgList(n_signal_total, n_back));
       
    RooPlot *frame = InvariantMass.frame(RooFit::Title("Invariant Mass"));
    RooPlot *frame_new = InvariantMass.frame(RooFit::Title("Invariant Mass"));
       
    double* output = new double[4];
    
    RooFitResult* fitres = new RooFitResult; //saves fit result
    fitres = model->fitTo(*dh, RooFit::Save());
    
    RooRealVar* yield_ALL = (RooRealVar*) fitres->floatParsFinal().find("n_signal_total");
    output[0] = yield_ALL->getVal();
    output[2] = yield_ALL->getError();
    
    c_all->cd();
    
    frame->SetTitle("ALL");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    //frame->SetYTitle(Form("Events / %3.1f MeV/c^{2}",Data_ALL->GetBinWidth(1)*1000));
    dh->plotOn(frame);
    model->plotOn(frame);
    frame->Draw("");
    // OUTPUT ARRAY
    
       
    model->plotOn(frame,RooFit::Components("GS"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame,RooFit::Components("CB"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    //model->paramOn(frame,RooFit::Layout(0.55, 0.99, 0.8));
    frame->Draw("");
    c_all->SaveAs(strcat(strcat(file, condition),all_pdf));
       
    frame_new->SetTitle("PASS");
    frame_new->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    //frame_new->SetYTitle(Form("Events / %3.1f MeV/c^{2}",Data_ALL->GetBinWidth(1)*1000));
    dh_pass->plotOn(frame_new);

    //PASSING PROBE CANVAS
    fitres = model->fitTo(*dh_pass, RooFit::Save());
    
    RooRealVar* yield_PASS = (RooRealVar*) fitres->floatParsFinal().find("n_signal_total");
    output[1] = yield_PASS->getVal();
    output[3] = yield_PASS->getError();
       
    c_pass->cd();

    dh_pass->plotOn(frame_new);
    model->plotOn(frame_new);
           
    model->plotOn(frame_new,RooFit::Components("GS"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame_new,RooFit::Components("CB"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame_new,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    //model->paramOn(frame_new,RooFit::Layout(0.55, 0.99, 0.8));
    frame_new->Draw("");
    c_pass->SaveAs(strcat(strcat(file,condition), pass_pdf));
    
    // DELETING ALLOCATED MEMORY
    delete file0;
    //delete DataTree; - Deleting TTree deletes TFile
    
    delete Data_ALL;
    delete Data_PASSING;
    //
    delete dh;
    delete dh_pass;
    //
    delete cutvar;
    delete redeuce;
    //
    delete signal;
    //
    delete c_all;
    delete c_pass;
    //
    delete model;
    //delete fitres;
    
    return output;
}

string* get_conditions(int bin_n, double* bins, string quantity = "ProbeMuon_Pt")
{
    string* conditions = new string[bin_n];
    for (int i = 0; i < bin_n; i++)
    {
        conditions[i] = quantity + ">" + to_string(bins[i]) + " && " + quantity + "<" + to_string(bins[i+1]);
    }
    return conditions;
}

TH1F* make_hist(string name, double** values, int qnt, int bin_n, Double_t* binning, bool IsDataMc, bool DRAW = false)
{
    //AddBinContent
    //HISTOGRAM NEEDS TO HAVE VARIABLE BINS
   
    TH1F* hist = new TH1F(name.c_str(), name.c_str(), bin_n, binning);

    for (int i = 0; i < bin_n; i++)
    {
        hist->SetBinContent(i, values[i][qnt]);
        if (IsDataMc == false)
            hist->SetBinError(i, values[i][qnt+2]);
    }
    if (DRAW)
    {
        TCanvas* xperiment = new TCanvas;
        xperiment->cd();
        hist->Draw();
    }
    return hist;
}

TEfficiency* get_efficiency(TH1F* ALL, TH1F* PASS)
{
    TFile* pFile = new TFile("Efficiency.root","recreate");
    TEfficiency* pEff = new TEfficiency();
    pEff->SetName("Efficiency");
    pEff->SetPassedHistogram(*PASS, "f");
    pEff->SetTotalHistogram (*ALL,"f");
    
    pEff->SetDirectory(gDirectory);
    pFile->Write();
    
    TCanvas* oi = new TCanvas();
    oi->cd();
    pEff->Draw();
    
    gPad->Update();

    //Set range in y axis

    auto graph = pEff->GetPaintedGraph();
    graph->SetMinimum(0.8);
    graph->SetMaximum(1.2);
    gPad->Update();
    
    return pEff;
}

double* McYield(string condition)
{
    TFile *file0    = TFile::Open("DATA/JPsi/T&P_JPSI_DATA_MC.root");
    TTree *DataTree = (TTree*)file0->Get(("JPSI_DATA"));
    
    //tenho de variar os valores deste corte para diferentes cortes em Pt
    double _mmin = 2.95;  double _mmax = 3.23;
    //double _mmin = 8.5;  double _mmax = 11;
    
    RooRealVar PassingProbeStandAloneMuon("PassingProbeStandAloneMuon", "PassingProbeStandAloneMuon", 0, 1);
    
    RooRealVar InvariantMass("InvariantMass", "InvariantMass", _mmin, _mmax);
    RooRealVar ProbeMuon_Pt("ProbeMuon_Pt", "ProbeMuon_Pt", 0, 60);
    RooRealVar ProbeMuon_Eta("ProbeMuon_Eta", "ProbeMuon_Eta", -3, 3);
    RooRealVar ProbeMuon_Phi("ProbeMuon_Phi", "ProbeMuon_Phi", -3.5, 3.5);
    
    RooFormulaVar* redeuce = new RooFormulaVar("PPTM", condition.c_str(), RooArgList(ProbeMuon_Pt));
    RooDataSet *Data_ALL    = new RooDataSet("DATA", "DATA", DataTree, RooArgSet(InvariantMass, PassingProbeStandAloneMuon, ProbeMuon_Pt),*redeuce);
    RooFormulaVar* cutvar = new RooFormulaVar("PPTM", strcat(condition.c_str(), "&& PassingProbeStandAloneMuon == 1"), RooArgList(PassingProbeStandAloneMuon, ProbeMuon_Pt));
    RooDataSet *Data_PASSING = new RooDataSet("DATA_PASS", "DATA_PASS", DataTree, RooArgSet(InvariantMass, PassingProbeStandAloneMuon, ProbeMuon_Pt), *cutvar);//
    
    double* output = new double[2];
    output[0] = Data_ALL->sumEntries();
    output[1] = Data_PASSING->sumEntries();
    return output;
}

//PROCURAR TEFFICIENCY
void JPsi_Fit()
{
    bool DataIsMC = false;
    
    double bins[] = {2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12, 12.5, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 40, 60};
    int bin_n = 105;
    //double bins[] = {-2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.5, 0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
    //int bin_n = 19;       //-- BINS USED TO CALCULATE PHI
    
    string* conditions = get_conditions(bin_n, bins);
    double ** yields_n_errs = new double*[bin_n]; // [yield_all, yield_pass, err_all, err_pass]
    
    for (int i = 0; i < bin_n; i++)
    {
        if (DataIsMC)
            yields_n_errs[i] = McYield(conditions[i]);
        else
            yields_n_errs[i] = doFit(conditions[i]);
    }
    
    
    TH1F *yield_ALL  = make_hist("ALL", yields_n_errs, 0, bin_n, bins, DataIsMC);
    TH1F *yield_PASS = make_hist("PASS", yields_n_errs, 1,bin_n, bins, DataIsMC);
    
     // saves histograms to .root
     TFile* EfficiencyFile = TFile::Open("Histograms.root","RECREATE");
     yield_ALL->SetDirectory(gDirectory);
     yield_PASS->SetDirectory(gDirectory);
     EfficiencyFile->Write();
     
     get_efficiency(yield_ALL, yield_PASS);
     //compare_efficiency();
     
     //delete bins;
     delete[] yields_n_errs;
     delete[] conditions;
     
     
}
