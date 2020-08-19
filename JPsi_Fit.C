#include "RooFitResult.h"

//EFFICIENCY TEST

//PROCURAR TEFFICIENCY
void JPsi_Fit()
{
    //ACCESSES THE DIRECTORY IN WHICH THE HISTOGRAM IS STORED
    TFile* files            = new TFile("DATA/JPsi_genData.root"); //PATH TO HISTOGRAM
    TDirectory* hist_files  = (TDirectory*)files->Get("histograms;1");
    
    
    //CREATES HISTOGRAMS
    TH1* hist_all  = (TH1*)hist_files->Get("AllMuonInvariantMass;1");
    TH1* hist_pass = (TH1*)hist_files->Get("PassingMuonInvariantMass;1");
    hist_all->Draw(); // CHECKS PLOT
 
    
    //FITTING PARAMETERS
    double _mmin = 2.85;  double _mmax = 3.3;
    
    
    double mass_peak1 = 3.1;
    
    TCanvas* c_all  = new TCanvas;
    TCanvas* c_pass = new TCanvas;
    
    //DECLARE OBSERVABLE X - INVARIANT MASS BETWEEN _mmin AND _mmax
    RooRealVar mass_all("mass_all","mass_all",_mmin,_mmax);
    
    // Create a binned dataset that imports contents of TH1 and associates itscontents to observable 'mass'
    RooDataHist dh("dh","dh",mass_all,RooFit::Import(*hist_all));
    RooDataHist dh_pass("dh_pass","dh_pass",mass_all,RooFit::Import(*hist_pass));


    
    RooRealVar lambda("lambda","lambda",-1.3,-10.,10.);
    RooExponential background("background", "background", mass_all, lambda);
    
    RooRealVar mean("mean","mean",_mmin,3.2);
    RooRealVar sigma("sigma","sigma",0.05*(_mmax-_mmin),0.,0.5*(_mmax-_mmin));
    RooRealVar alpha("alpha", "alpha", 1., 0.01, 3.);
    RooRealVar n("n", "n", 1.00);
    n.setConstant(kTRUE);
    
    //FIT FUNCTIONS
    RooGaussian gaussian("GS","GS",mass_all,mean,sigma);
    RooCBShape crystalball("CB", "CB", mass_all, mean, sigma, alpha, n);
    
    
    double n_signal_initial_total =(dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak1)) -dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>.015",mass_peak1,mass_peak1))) / dh.sumEntries();
    
    //signal = gaussian1*frac1 + gaussian2*frac2 + gaussian3*(1-(frac1 + frac2))
    //S(signal)d mass = 1
    RooRealVar frac("frac","frac",0.5,0.,1.);
    //para RooArgList N-1, assume como frações
    RooAddPdf* signal;
    signal = new RooAddPdf("signal", "signal", RooArgList(gaussian, crystalball), RooArgList(frac));
    
    double n_back_initial = 1. - n_signal_initial_total;
    
    RooRealVar n_signal_total("n_signal_total","n_signal_total",n_signal_initial_total,0.,dh.sumEntries());
    
    RooRealVar n_back("n_back","n_back",n_back_initial,0.,dh.sumEntries());
    //para RooArgList N, assume como normalizações
    //modelo_total = n_signal_total*signal + n_back*background
    RooAddPdf* model;
    model = new RooAddPdf("model","model", RooArgList(*signal, background),RooArgList(n_signal_total, n_back));
    
    RooPlot *frame = mass_all.frame(RooFit::Title("Invariant Mass"));
    RooPlot *frame_new = mass_all.frame(RooFit::Title("Invariant Mass"));
    
    RooFitResult* fitres = new RooFitResult; //saves fit result
    fitres = model->fitTo(dh);
    
    //PLOTTING
    
    //legends
    
    //FAKE FUNCTIONS FOR LEGEND CREATION
    TF1* fake1 = new TF1();
    TF1* fake2 = new TF1();
    TF1* fake3 = new TF1();
    TF1* fake5 = new TF1();
    
    fake1->SetLineColor(kBlue);
    fake2->SetLineColor(kGreen);
    fake3->SetLineColor(kMagenta - 5);
    fake5->SetLineColor(kRed);
    
    fake2->SetLineStyle(kDashDotted);
    fake3->SetLineStyle(kDashDotted);
    fake5->SetLineStyle(kDashDotted);
    
    TLegend* tl = new TLegend(0.70,0.70,0.90,0.92);
    tl->SetTextSize(0.04);
    tl->AddEntry(model,"Data", "p");
    tl->AddEntry(fake1,"Total Fit", "l");
    tl->AddEntry(fake2,"Gaussian", "l");
    tl->AddEntry(fake3,"CrystalBall", "l");
    tl->AddEntry(fake5,"background", "l");
    
    //ALL EVENTS CANVAS
    c_all->cd();
        
    frame->SetTitle("ALL");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    frame->SetYTitle(Form("Events / %3.1f MeV/c^{2}",hist_all->GetBinWidth(1)*1000));
    dh.plotOn(frame);
            
    model->plotOn(frame);
            
    model->plotOn(frame,RooFit::Components("GS"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame,RooFit::Components("CB"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    //model->paramOn(frame,RooFit::Layout(0.55, 0.99, 0.8));
    frame->Draw("");
    tl->Draw();
    c_all->SaveAs("Result/invariant_mass_ALL.pdf");
    
    frame_new->SetTitle("PASS");
    frame_new->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    frame_new->SetYTitle(Form("Events / %3.1f MeV/c^{2}",hist_all->GetBinWidth(1)*1000));
    dh_pass.plotOn(frame_new);

    //PASSING PROBE CANVAS
    model->fitTo(dh_pass);
    
    c_pass->cd();

    dh_pass.plotOn(frame_new);
    model->plotOn(frame_new);
        
    model->plotOn(frame_new,RooFit::Components("GS"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame_new,RooFit::Components("CB"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame_new,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    //model->paramOn(frame_new,RooFit::Layout(0.55, 0.99, 0.8));
    frame_new->Draw("");
    tl->Draw();
    c_pass->SaveAs("Result/invariant_mass_PASS.pdf");
}
