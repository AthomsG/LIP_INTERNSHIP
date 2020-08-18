#include "RooFitResult.h"

//EFFICIENCY TEST

//PROCURAR TEFFICIENCY
void JPsi_fit()
{
    //ACCESSES THE DIRECTORY IN WHICH THE HISTOGRAM IS STORED
    TFile* files            = new TFile("generated_hist.root"); //PATH TO HISTOGRAM
    TDirectory* hist_files  = (TDirectory*)files->Get("histograms;1");
    
    //CREATES HISTOGRAMS
    TH1* hist_all  = (TH1*)hist_files->Get("AllMuonInvariantMass;1");
    TH1* hist_pass = (TH1*)hist_files->Get("Passing trackerMuonInvariantMass;1");
    hist_all->Draw(); // CHECKS PLOT
 
    //FITTING PARAMETERS
    double _mmin = 6.1;  double _mmax = 6.4;
    
    
    double mass_peak1 = 9.46030;
    double mass_peak2 = 10.02326;
    double mass_peak3 = 10.3552;
    
    TCanvas* c_all  = new TCanvas;
    TCanvas* c_pass = new TCanvas;
    
    //DECLARE OBSERVABLE X - INVARIANT MASS BETWEEN _mmin AND _mmax
    RooRealVar mass_all("mass_all","mass_all",_mmin,_mmax);
    
    // Create a binned dataset that imports contents of TH1 and associates itscontents to observable 'mass'
    RooDataHist dh("dh","dh",mass_all,RooFit::Import(*hist_all));
    RooDataHist dh_pass("dh_pass","dh_pass",mass_all,RooFit::Import(*hist_pass));


    
    RooRealVar lambda("lambda","lambda",-1.3,-10.,10.);
    RooExponential background("background", "background", mass_all, lambda);
    
    RooRealVar sigma("sigma","sigma",0.05*(_mmax-_mmin),0.,0.5*(_mmax-_mmin));
    
    RooRealVar mean1("mean1","mean1",_mmin,(mass_peak1+mass_peak2)/2.);
    RooRealVar mean2("mean2","mean2",(mass_peak1+mass_peak2)/2.,(mass_peak3+mass_peak2)/2.);
    RooRealVar mean3("mean3","mean3",(mass_peak3+mass_peak2)/2.,_mmax);
    
    //FIT FUNCTIONS
    
    // --Gaussian as the signal pdf
    RooGaussian gaussian1("signal1","signal1",mass_all,mean1,sigma);
    RooGaussian gaussian2("signal2","signal2",mass_all,mean2,sigma);
    RooGaussian gaussian3("signal3","signal3",mass_all,mean3,sigma);
    
    double n_signal_initial1 =(dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak1)) -dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>.015",mass_peak1,mass_peak1))) / dh.sumEntries();
    double n_signal_initial2 =(dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak2)) -dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>.015",mass_peak2,mass_peak2))) / dh.sumEntries();
    double n_signal_initial3 =(dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak3)) -dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>.015",mass_peak3,mass_peak3))) / dh.sumEntries();
    
    double n_signal_initial_total = n_signal_initial1 + n_signal_initial2 + n_signal_initial3;
    
    //signal = gaussian1*frac1 + gaussian2*frac2 + gaussian3*(1-(frac1 + frac2))
    //S(signal)d mass = 1
    RooRealVar frac1("frac1","frac1",0.333,0.,1.);
    RooRealVar frac2("frac2","frac2",0.333,0.,1.);
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
    TF1* fake4 = new TF1();
    TF1* fake5 = new TF1();
    
    fake1->SetLineColor(kBlue);
    fake2->SetLineColor(kGreen);
    fake3->SetLineColor(kMagenta - 5);
    fake4->SetLineColor(kOrange);
    fake5->SetLineColor(kRed);
    
    fake2->SetLineStyle(kDashDotted);
    fake3->SetLineStyle(kDashDotted);
    fake4->SetLineStyle(kDashDotted);
    fake5->SetLineStyle(kDashDotted);
    
    TLegend* tl = new TLegend(0.70,0.70,0.90,0.92);
    tl->SetTextSize(0.04);
    tl->AddEntry(model,"Data", "p");
    tl->AddEntry(fake1,"Total Fit", "l");
    tl->AddEntry(fake2,"Peak 1", "l");
    tl->AddEntry(fake3,"Peak 2", "l");
    tl->AddEntry(fake4,"Peak 3", "l");
    tl->AddEntry(fake5,"background", "l");
    
    //ALL EVENTS CANVAS
    c_all->cd();
        
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
        
    model->plotOn(frame_new,RooFit::Components("signal1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame_new,RooFit::Components("signal2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame_new,RooFit::Components("signal3"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
    model->plotOn(frame_new,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    //model->paramOn(frame_new,RooFit::Layout(0.55, 0.99, 0.8));
    frame_new->Draw("");
    tl->Draw();
    c_pass->SaveAs("Result/invariant_mass_PASS.pdf");
}
