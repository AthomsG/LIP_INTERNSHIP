using namespace RooFit ;

void Fit_From_Run()
{
    TFile *file0  = TFile::Open("Run2011_upsilon.root");
    TTree *DataTree = (TTree*)file0->Get(("tagandprobe/AnalysisTree"));
    
    double _mmin = 8.7;  double _mmax = 11;

    RooRealVar InvariantMass("InvariantMass", "InvariantMass", _mmin, _mmax);
    
    RooRealVar PassingProbeTrackingMuon("PassingProbeTrackingMuon", "PassingProbeTrackingMuon", 0, 1);
    RooRealVar PassingProbeStandAloneMuon("PassingProbeStandAloneMuon", "PassingProbeStandAloneMuon", 0, 1);
    RooRealVar PassingProbeGlobalMuon("PassingProbeGlobalMuon", "PassingProbeGlobalMuon", 0, 1);
    
    RooDataSet *Data_ALL = new RooDataSet("DATA", "DATA", DataTree, RooArgSet(InvariantMass)/*, "InvariantMass > 9.5"*/);
    
    // MUDAR SIGNAL PARA CRYSTALLBALL
    // APLICAR CORTE A GRANDEZA PEDIDA (i.e "pt<5 and pt>5") e retirar o yield, para bins
    //RooDataSet* Data = Data_ALL.Reduce("");
    
    RooFormulaVar* condition = new RooFormulaVar("condition", "PassingProbeTrackingMuon == 1", RooArgList(PassingProbeTrackingMuon, PassingProbeStandAloneMuon, PassingProbeGlobalMuon));
    //RooDataSet *Data_pass = new RooDataSet("DATA_pass", "DATA_pass", DataTree, RooArgSet(InvariantMass), condition);
    
    double mass_peak1 = 9.46030;
    double mass_peak2 = 10.02326;
    double mass_peak3 = 10.3552;
    
    TCanvas* c_all  = new TCanvas;
    RooPlot *frame = InvariantMass.frame(RooFit::Title("Invariant Mass"));
    
    // Create an unbinned dataset that imports contents of TH1 and associates itscontents to observable 'mass'
    RooDataHist* dh = Data_ALL->binnedClone("All Invariant Mass", "All Invariant Mass");
    
    // PROBLEMAS NO BACKGROUND
    RooRealVar a0("a0", "a0", 0., -100, 100);
    RooRealVar a1("a1", "a1", 0., -100, 100);
    //RooRealVar a2("a2", "a2", -1000, 1000);
    RooChebychev background("cpol","cpol", InvariantMass, RooArgList(a0,a1));
    
    //RooRealVar lambda("lambda","lambda",-1.3,-10.,10.);
    //RooExponential background("background", "background", InvariantMass, lambda);
    
    RooRealVar sigma("sigma","sigma",0.08,0.05,0.1);
    RooRealVar mean1("mean1","mean1",9.4,9.3,9.5);
    RooRealVar mean2("mean2","mean2",10, 9.9, 10.1);
    RooRealVar mean3("mean3","mean3",10.35, 10.2, 10.4);
    //FIT FUNCTIONS
    
    // --Gaussian as the signal pdf
    RooGaussian gaussian1("signal1","signal1",InvariantMass,mean1,sigma);
    RooGaussian gaussian2("signal2","signal2",InvariantMass,mean2,sigma);
    RooGaussian gaussian3("signal3","signal3",InvariantMass,mean3,sigma);
    
    double n_signal_initial1 =(dh->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak1)) -dh->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak1,mass_peak1))) / dh->sumEntries();
    double n_signal_initial2 =(dh->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak2)) -dh->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak2,mass_peak2))) / dh->sumEntries();
    double n_signal_initial3 =(dh->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak3)) -dh->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak3,mass_peak3))) / dh->sumEntries();
    
    double n_signal_initial_total = n_signal_initial1 + n_signal_initial2 + n_signal_initial3;
    
    //signal = gaussian1*frac1 + gaussian2*frac2 + gaussian3*(1-(frac1 + frac2))
    //S(signal)d mass = 1
    RooRealVar frac1("frac1","frac1",0.333,-1.,1.);
    RooRealVar frac2("frac2","frac2",0.333,-1.,1.);
    //para RooArgList N-1, assume como frações
    RooAddPdf* signal;
    signal = new RooAddPdf("signal", "signal", RooArgList(gaussian1, gaussian2,gaussian3), RooArgList(frac1, frac2));
    
    double n_back_initial = 1. - n_signal_initial1 - n_signal_initial2 -n_signal_initial3;
    
    RooRealVar n_signal_total("n_signal_total","n_signal_total",n_signal_initial_total,0.,dh->sumEntries());
    
    RooRealVar n_back("n_back","n_back",n_back_initial,0.,dh->sumEntries());
    //para RooArgList N, assume como normalizações
    //modelo_total = n_signal_total*signal + n_back*background
    RooAddPdf* model;
    model = new RooAddPdf("model","model", RooArgList(*signal, background),RooArgList(n_signal_total, n_back));
    
    RooFitResult* fitres = new RooFitResult; //saves fit result
    fitres = model->fitTo(*dh);
    
    fitres = model->fitTo(*Data_ALL);
    
    //RooRealVar* pass_mean1 = (RooRealVar*) fitres->floatParsFinal().find("mean1");
    //RooRealVar* pass_mean2 = (RooRealVar*) fitres->floatParsFinal().find("mean2");
    //RooRealVar* pass_mean3 = (RooRealVar*) fitres->floatParsFinal().find("mean3");
    //
    //RooRealVar* pass_sigma = (RooRealVar*) fitres->floatParsFinal().find("sigma");
    //
    //double mean1_value = pass_mean1->getVal();
    //double mean2_value = pass_mean2->getVal();
    //double mean3_value = pass_mean3->getVal();
    //
    //double sigma_value = pass_sigma->getVal();
    
    frame->SetTitle("ALL");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    //frame->SetYTitle(Form("Events / %3.1f MeV/c^{2}",dh->GetBinWidth(1)*1000));
    dh->plotOn(frame);
    
    
    model->plotOn(frame);
    model->plotOn(frame,RooFit::Components("signal1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame,RooFit::Components("signal2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame,RooFit::Components("signal3"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
    model->plotOn(frame,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    //model->paramOn(frame,RooFit::Layout(0.55, 0.99, 0.8));
    
    frame->Draw("");
    
    c_all->SaveAs("fit_DATA.pdf");
}
