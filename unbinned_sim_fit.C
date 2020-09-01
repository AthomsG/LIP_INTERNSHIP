using namespace RooFit ;

void unbinned_sim_fit()
{
    TFile *file0  = TFile::Open("DATA/T&P_UPSILON_DATA.root");
    TTree *DataTree = (TTree*)file0->Get(("UPSILON_DATA"));
    
    //I need to add quantities branches to DataTree
    double _mmin = 8.7;  double _mmax = 11;

    RooRealVar PassingProbeTrackingMuon("PassingProbeTrackingMuon", "PassingProbeTrackingMuon", 0, 1);
    
    RooRealVar InvariantMass("InvariantMass", "InvariantMass", _mmin, _mmax);
    RooRealVar ProbeMuon_Pt("ProbeMuon_Pt", "ProbeMuon_Pt", 0, 37.2);
    RooRealVar ProbeMuon_Eta("ProbeMuon_Eta", "ProbeMuon_Eta", -3, 3);
    RooRealVar ProbeMuon_Phi("ProbeMuon_Phi", "ProbeMuon_Phi", -3.5, 3.5);
    
    RooFormulaVar* redeuce = new RooFormulaVar("PPTM", "ProbeMuon_Pt < 4 && ProbeMuon_Pt > 0", RooArgList(ProbeMuon_Pt));
    RooDataSet *Data_ALL    = new RooDataSet("DATA", "DATA", DataTree, RooArgSet(InvariantMass, PassingProbeTrackingMuon, ProbeMuon_Pt), *redeuce);
    //RooFormulaVar* cutvar = new RooFormulaVar("PPTM", "PassingProbeTrackingMuon ==  1", RooArgList(PassingProbeTrackingMuon));
    RooFormulaVar* cutvar = new RooFormulaVar("PPTM", "ProbeMuon_Pt < 4 && ProbeMuon_Pt > 0 && PassingProbeTrackingMuon == 1", RooArgList(ProbeMuon_Pt, PassingProbeTrackingMuon));
    RooDataSet *Data_PASSING = new RooDataSet("DATA_PASS", "DATA_PASS", DataTree, RooArgSet(InvariantMass, PassingProbeTrackingMuon, ProbeMuon_Pt), *cutvar);//
    
    //BINNING DATASET
    RooDataHist* dh_ALL     = Data_ALL->binnedClone();
    RooDataHist* dh_PASSING = Data_PASSING->binnedClone();
    
    // APLICAR CORTE A GRANDEZA PEDIDA (i.e "pt<5 and pt>5") e retirar o yield, para bins
    //RooDataSet* Data = Data_ALL.Reduce("");
    
    double mass_peak1 = 9.46030;
    double mass_peak2 = 10.02326;
    double mass_peak3 = 10.3552;
    
    TCanvas* c_all  = new TCanvas;
    TCanvas* c_pass = new TCanvas;
    
    RooPlot *frame = InvariantMass.frame(RooFit::Title("Invariant Mass"));

    // BACKGROUND VARIABLES
    RooRealVar a0("a0", "a0", 0, -100, 100);
    RooRealVar a1("a1", "a1", 0, -100, 100);
    
    RooRealVar a0_pass("a0_pass", "a0_pass", 0., -100, 100);
    RooRealVar a1_pass("a1_pass", "a1_pass", 0., -100, 100);
    // BACKGROUND FUNCTION
    RooChebychev background("background","background", InvariantMass, RooArgList(a0,a1));
    
    RooChebychev background_pass("cpol_pass","cpol_pass", InvariantMass, RooArgList(a0_pass,a1_pass));
    // GAUSSIAN VARIABLES
    RooRealVar sigma("sigma","sigma",0.08,0.05,0.1);
    RooRealVar mean1("mean1","mean1",mass_peak1,9.4,9.47);
    RooRealVar mean2("mean2","mean2",mass_peak2, 10, 10.03);
    RooRealVar mean3("mean3","mean3",mass_peak3, 10.3, 10.4);
    // CRYSTAL BALL VARIABLES
    RooRealVar alpha("alpha","alpha", 1.4384e+00, 1.43, 1.44);
    RooRealVar n("n", "n", 1.6474e+01, 16., 17.);
    // FIT FUNCTIONS
    RooCBShape  gaussian1("signal1","signal1",InvariantMass,mean1,sigma, alpha, n);
    RooGaussian gaussian2("signal2","signal2",InvariantMass,mean2,sigma);
    RooGaussian gaussian3("signal3","signal3",InvariantMass,mean3,sigma);
    
    double n_signal_initial1 =(Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak1)) -Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak1,mass_peak1))) / Data_ALL->sumEntries();
    double n_signal_initial2 =(Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak2)) -Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak2,mass_peak2))) / Data_ALL->sumEntries();
    double n_signal_initial3 =(Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak3)) -Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak3,mass_peak3))) / Data_ALL->sumEntries();
    
    double n_signal_initial_total = n_signal_initial1 + n_signal_initial2 + n_signal_initial3;
    
    RooRealVar frac1("frac1","frac1",7.1345e-01,0.7134,0.7135);
    RooRealVar frac2("frac2","frac2",1.9309e-01,0.193,0.194);
 
    RooAddPdf* signal;
    RooAddPdf* signal_pass;
    
    signal      = new RooAddPdf("signal", "signal", RooArgList(gaussian1, gaussian2,gaussian3), RooArgList(frac1, frac2));

    double n_back_initial = 1. - n_signal_initial1 - n_signal_initial2 -n_signal_initial3;
    
    RooRealVar n_signal_total("n_signal_total","n_signal_total",n_signal_initial_total,0.,Data_ALL->sumEntries());
    RooRealVar n_signal_total_pass("n_signal_total_pass","n_signal_total_pass",n_signal_initial_total,0.,Data_PASSING->sumEntries());
    
    RooRealVar n_back("n_back","n_back",n_back_initial,0.,Data_ALL->sumEntries());
    RooRealVar n_back_pass("n_back_pass","n_back_pass",n_back_initial,0.,Data_PASSING->sumEntries());

    RooAddPdf* model;
    RooAddPdf* model_pass;
    
    model      = new RooAddPdf("model","model", RooArgList(*signal, background),RooArgList(n_signal_total, n_back));
    model_pass = new RooAddPdf("model_pass", "model_pass", RooArgList(*signal, background),RooArgList(n_signal_total_pass, n_back_pass));
    
    // SIMULTANEOUS FIT
    
    // C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s
    // ---------------------------------------------------------------------------
    RooCategory sample("sample","sample") ;
    sample.defineType("All") ;
    sample.defineType("PASSING") ;
    
    //WHEN DOING UNBINNED, CHANGE TO DATASET
    RooDataHist combData("combData","combined data",InvariantMass,Index(sample),Import("ALL",*dh_ALL),Import("PASSING",*dh_PASSING));
   
    // C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
    // -----------------------------------------------------------------------------------
    
    RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
   
    simPdf.addPdf(*model,"ALL");
    simPdf.addPdf(*model_pass,"PASSING");
    
    // TAKES REALLY LONG TO FIT THIS SHIT
    RooFitResult* fitres = new RooFitResult;
    fitres = simPdf.fitTo(combData, RooFit::Save());
    
    fitres->Print();
         
    // Some Random Stuff like... whaaaa....
    frame->SetTitle("ALL");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    //frame->SetYTitle(Form("Events / %3.1f MeV/c^{2}",dh_ALL->GetBinWidth(1)*1000));
    Data_ALL->plotOn(frame);
    
    model->plotOn(frame);
    model->plotOn(frame,RooFit::Components("signal1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame,RooFit::Components("signal2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame,RooFit::Components("signal3"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
    model->plotOn(frame,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    //model->paramOn(frame,RooFit::Layout(0.55, 0.99, 0.8));
    
    c_all->cd();
    frame->Draw("");
   
    //TESTING THE NEW DATASET
    RooPlot *frame_pass = InvariantMass.frame(RooFit::Title("Invariant Mass"));
    
    c_pass->cd();
    
    frame_pass->SetTitle("PASSING");
    frame_pass->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    Data_PASSING->plotOn(frame_pass);
    
    model_pass->plotOn(frame_pass);
    model_pass->plotOn(frame_pass,RooFit::Components("signal1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model_pass->plotOn(frame_pass,RooFit::Components("signal2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model_pass->plotOn(frame_pass,RooFit::Components("signal3"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
    model_pass->plotOn(frame_pass,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    
    frame_pass->Draw();
    
    c_pass->SaveAs("Result/fit_DATA_PASS.pdf");
    c_all->SaveAs("Result/fit_DATA.pdf");
}

