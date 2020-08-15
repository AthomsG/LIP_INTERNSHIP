//EFFICIENCY TEST
{
    //ACCESSES THE DIRECTORY IN WHICH THE HISTOGRAM IS STORED
    TFile* files            = new TFile("generated_hist.root"); //PATH TO HISTOGRAM
    TDirectory* hist_files  = (TDirectory*)files->Get("histograms;1");
    
    //CREATES HISTOGRAMS
    TH1* hist_all  = (TH1*)hist_files->Get("AllMuonInvariantMass;1");
    //TH1* hist_all = (TH1*)hist_files->Get("Passing trackerMuonInvariantMass;1");
    //hist->Draw(); // CHECKS PLOT
 
    //FITTING PARAMETERS
    Double_t mmin(9), mmax(11);
    _mmin = mmin;  _mmax = mmax;
    
    double mass_peak1 = 9.46030;
    double mass_peak2 = 10.02326;
    double mass_peak3 = 10.3552;
    
    // I WASN'T ABLE TO USE THE SAME MODEL ON BOTH HISTOGRAMS SO I HAD TO DUPLICATE THIS CODE TWICE------------------------
    
        //DECLARE OBSERVABLE X - INVARIANT MASS BETWEEN _mmin AND _mmax
        RooRealVar mass_all("mass_all","mass_all",_mmin,_mmax);
        
        // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'mass'
        RooDataHist dh("dh","dh",mass_all,RooFit::Import(*hist_all));
        
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
        
        double n_signal_initial1 = (dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak1)) - dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>0.015",mass_peak1,mass_peak1))) / dh.sumEntries();
        double n_signal_initial2 = (dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak2)) - dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>0.015",mass_peak2,mass_peak2))) / dh.sumEntries();
        double n_signal_initial3 = (dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak3)) - dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>0.015",mass_peak3,mass_peak3))) / dh.sumEntries();
        
        double n_signal_initial_total = n_signal_initial1 + n_signal_initial2 + n_signal_initial3;
        
        RooRealVar frac1("frac1","frac1",0.333,0.,1.);
        RooRealVar frac2("frac2","frac2",0.333,0.,1.);

        RooAddPdf* signal;
        signal = new RooAddPdf("signal", "signal", RooArgList(gaussian1, gaussian2, gaussian3), RooArgList(frac1, frac2));


        double n_back_initial = 1. - n_signal_initial1 - n_signal_initial2 - n_signal_initial3;

        RooRealVar n_signal1("n_signal1","n_signal1",n_signal_initial1,0.,dh.sumEntries());
        RooRealVar n_signal2("n_signal2","n_signal2",n_signal_initial2,0.,dh.sumEntries());
        RooRealVar n_signal3("n_signal3","n_signal3",n_signal_initial3,0.,dh.sumEntries());
        RooRealVar n_signal_total("n_signal_total","n_signal_total",n_signal_initial_total,0.,dh.sumEntries());
        

    RooRealVar n_back("n_back","n_back",n_back_initial,0.,dh.sumEntries());


        RooAddPdf* model;
        model = new RooAddPdf("model","model", RooArgList(*signal, background), RooArgList(n_signal_total, n_back));

        RooPlot *frame = mass_all.frame(RooFit::Title("Invariant Mass"));
        
        model->fitTo(dh);
        
        dh.plotOn(frame);
        
        model->plotOn(frame);
        
        model->plotOn(frame,RooFit::Components("signal1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
        model->plotOn(frame,RooFit::Components("signal2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
        model->plotOn(frame,RooFit::Components("signal3"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
        model->plotOn(frame,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
        //model->paramOn(frame,RooFit::Layout(0.55, 0.99, 0.8));
    
    
    
    
    
    
// HERE IS THE REPLICATE----------------------------------------------------------------------------------------------------
    
    
    
    
    
        
         RooRealVar mass_pass("mass_pass","mass_pass",_mmin,_mmax);
           
           // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'mass'
           RooDataHist dh_pass("dh_pass","dh_pass",mass_pass,RooFit::Import(*hist_pass));
           
           RooRealVar lambda_pass("lambda_pass","lambda_pass",-1.3,-10.,10.);
           RooExponential background_pass("background_pass", "background_pass", mass_pass, lambda_pass);
           
           RooRealVar sigma_pass("sigma_pass","sigma_pass",0.05*(_mmax-_mmin),0.,0.5*(_mmax-_mmin));
           
           RooRealVar mean1_pass("mean1_pass","mean1_pass",_mmin,(mass_peak1+mass_peak2)/2.);
           RooRealVar mean2_pass("mean2_pass","mean2_pass",(mass_peak1+mass_peak2)/2.,(mass_peak3+mass_peak2)/2.);
           RooRealVar mean3_pass("mean3_pass","mean3_pass",(mass_peak3+mass_peak2)/2.,_mmax);
           
           //FIT FUNCTIONS
           
           // -- gaussian as the signal pdf
           RooGaussian gaussian_pass1("signal1_pass","signal1_pass",mass_pass,mean1_pass,sigma_pass);
           RooGaussian gaussian_pass2("signal2_pass","signal2_pass",mass_pass,mean2_pass,sigma_pass);
           RooGaussian gaussian_pass3("signal3_pass","signal3_pass",mass_pass,mean3_pass,sigma_pass);
           
           double n_signal_initial_pass1 = (dh_pass.sumEntries(TString::Format("abs(mass_pass-%g)<0.015",mass_peak1)) - dh_pass.sumEntries(TString::Format("abs(mass_pass-%g)<0.030&&abs(mass_pass-%g)>0.015",mass_peak1,mass_peak1))) / dh_pass.sumEntries();
           double n_signal_initial_pass2 = (dh_pass.sumEntries(TString::Format("abs(mass_pass-%g)<0.015",mass_peak2)) - dh_pass.sumEntries(TString::Format("abs(mass_pass-%g)<0.030&&abs(mass_pass-%g)>0.015",mass_peak2,mass_peak2))) / dh_pass.sumEntries();
           double n_signal_initial_pass3 = (dh_pass.sumEntries(TString::Format("abs(mass_pass-%g)<0.015",mass_peak3)) - dh_pass.sumEntries(TString::Format("abs(mass_pass-%g)<0.030&&abs(mass_pass-%g)>0.015",mass_peak3,mass_peak3))) / dh_pass.sumEntries();
           
           double n_signal_initial_pass_total = n_signal_initial_pass1 + n_signal_initial_pass2 + n_signal_initial_pass3;
           
           RooRealVar frac_pass1("frac_pass1","frac_pass1",0.333,0.,1.);
           RooRealVar frac_pass2("frac_pass2","frac_pass2",0.333,0.,1.);

           RooAddPdf* signal_pass;
           signal_pass = new RooAddPdf("signal_pass", "signal_pass", RooArgList(gaussian_pass1, gaussian_pass2, gaussian_pass3), RooArgList(frac_pass1, frac_pass2));


           double n_back_initial_pass = 1. - n_signal_initial_pass1 - n_signal_initial_pass2 - n_signal_initial_pass3;

           RooRealVar n_signal1_pass("n_signal1_pass","n_signal1_pass",n_signal_initial_pass1,0.,dh_pass.sumEntries());
           RooRealVar n_signal2_pass("n_signal2_pass","n_signal2_pass",n_signal_initial_pass2,0.,dh_pass.sumEntries());
           RooRealVar n_signal3_pass("n_signal3_pass","n_signal3_pass",n_signal_initial_pass3,0.,dh_pass.sumEntries());
           RooRealVar n_signal_total_pass("n_signal_total_pass","n_signal_total_pass",n_signal_initial_pass_total, 0.,dh_pass.sumEntries());
           //IS AT LOWER LIMLIT!? WHY IS IT FIXED ON 0?
           

           RooRealVar n_back_pass("n_back_pass","n_back_pass",n_back_initial_pass, -0.3,dh_pass.sumEntries());


           RooAddPdf* model_pass;
           model_pass = new RooAddPdf("model_pass","model_pass", RooArgList(*signal_pass, background_pass), RooArgList(n_signal_total_pass, n_back_pass));

           RooPlot *frame_pass = mass_pass.frame(RooFit::Title("Invariant Mass"));
           
           model_pass->fitTo(dh_pass);
           
           dh_pass.plotOn(frame_pass);
           
           model_pass->plotOn(frame_pass);
           
           model_pass->plotOn(frame_pass,RooFit::Components("signal1_pass"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
           model_pass->plotOn(frame_pass,RooFit::Components("signal2_pass"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
           model_pass->plotOn(frame_pass,RooFit::Components("signal3_pass"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
           model_pass->plotOn(frame_pass,RooFit::Components("background_pass"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    
    //FAKE FUNCTIONS FOR COLORS ---------------------------- AESTHETICS

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
    tl->AddEntry(model_pass, "Data", "p");
    tl->AddEntry(fake1,"Total Fit", "l");
    tl->AddEntry(fake2,"Peak 1", "l");
    tl->AddEntry(fake3,"Peak 2", "l");
    tl->AddEntry(fake4,"Peak 3", "l");
    tl->AddEntry(fake5,"background", "l");
    
    TCanvas* c_all  = new TCanvas;
    TCanvas* c_pass = new TCanvas;
    
    c_all->cd();
    
    frame->SetTitle("ALL");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    frame->SetYTitle(Form("Events / %3.1f MeV/c^{2}",hist_all->GetBinWidth(1)*1000));
    frame->Draw("");
    tl->Draw();
    
    c_pass->cd();
    
    frame_pass->SetTitle("PASSING");
    frame_pass->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    frame_pass->SetYTitle(Form("Events / %3.1f MeV/c^{2}",hist_all->GetBinWidth(1)*1000));
    frame_pass->Draw("");
    tl->Draw();
    
    //roofit_canvas.SaveAs("Slide.pdf");
}
