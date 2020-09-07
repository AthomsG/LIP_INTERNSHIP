using namespace RooFit ;

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

double* doFit(string condition, bool save = TRUE) // RETURNS ARRAY WITH [yield_all, yield_pass, err_all, err_pass]    -> OUTPUT ARRAY
{
    TFile *file0    = TFile::Open("DATA/T&P_UPSILON_DATA.root");
    TTree *DataTree = (TTree*)file0->Get(("UPSILON_DATA"));
    
    //tenho de variar os valores deste corte para diferentes cortes em Pt
    double _mmin = 9.1;  double _mmax = 10.6;
    //double _mmin = 8.5;  double _mmax = 11;
    
    RooRealVar PassingProbeTrackingMuon("PassingProbeTrackingMuon", "PassingProbeTrackingMuon", 0, 1);
    
    RooRealVar InvariantMass("InvariantMass", "InvariantMass", _mmin, _mmax);
    RooRealVar ProbeMuon_Pt("ProbeMuon_Pt", "ProbeMuon_Pt", 0, 60);
    RooRealVar ProbeMuon_Eta("ProbeMuon_Eta", "ProbeMuon_Eta", -3, 3);
    RooRealVar ProbeMuon_Phi("ProbeMuon_Phi", "ProbeMuon_Phi", -3.5, 3.5);
    
    RooFormulaVar* redeuce = new RooFormulaVar("PPTM", condition.c_str(), RooArgList(ProbeMuon_Pt));
    RooDataSet *Data_ALL    = new RooDataSet("DATA", "DATA", DataTree, RooArgSet(InvariantMass, PassingProbeTrackingMuon, ProbeMuon_Pt),*redeuce);
    RooFormulaVar* cutvar = new RooFormulaVar("PPTM", strcat(condition.c_str(), "&& PassingProbeTrackingMuon == 1"), RooArgList(PassingProbeTrackingMuon, ProbeMuon_Pt));
    RooDataSet *Data_PASSING = new RooDataSet("DATA_PASS", "DATA_PASS", DataTree, RooArgSet(InvariantMass, PassingProbeTrackingMuon, ProbeMuon_Pt), *cutvar);//
    
    //BINNING DATASET
    //DE RooDataSet -> TH1 -> RooDataHist
    
    RooDataHist* dh_ALL     = Data_ALL->binnedClone();
    RooDataHist* dh_PASSING = Data_PASSING->binnedClone();
    
    //TH1F* dh_ALL     = Data_ALL->CreateHistogram("dh_ALL", InvariantMass)
    //RooDataHist* dh_PASSING = Data_PASSING->binnedClone();
    
    // APLICAR CORTE A GRANDEZA PEDIDA (i.e "pt<5 and pt>5") e retirar o yield, para bins
    double mass_peak1 = 9.46030;
    double mass_peak2 = 10.02326;
    double mass_peak3 = 10.3552;
    
    TCanvas* c_all  = new TCanvas;
    TCanvas* c_pass = new TCanvas;
    
    RooPlot *frame = InvariantMass.frame(RooFit::Title("Invariant Mass"));
    // BACKGROUND VARIABLES
    RooRealVar a0("a0", "a0", 2.5875e-02, -1., 1.);
    RooRealVar a1("a1", "a1", -7.8407e-02, -1., 1.);

    // BACKGROUND FUNCTION
    RooChebychev background("background","background", InvariantMass, RooArgList(a0,a1));
    
    // GAUSSIAN VARIABLES
    RooRealVar sigma("sigma","sigma",0.08,0.05,0.11); //0.1
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
    
    double n_signal_initial1 =(Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak1))-Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak1,mass_peak1))) / Data_ALL->sumEntries();
    double n_signal_initial2 =(Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak2))-Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak2,mass_peak2))) / Data_ALL->sumEntries();
    double n_signal_initial3 =(Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.015",mass_peak3))-Data_ALL->sumEntries(TString::Format("abs(InvariantMass-%g)<0.030&&abs(InvariantMass-%g)>.015",mass_peak3,mass_peak3))) / Data_ALL->sumEntries();
    
    double n_signal_initial_total = n_signal_initial1 + n_signal_initial2 + n_signal_initial3;
    
    RooRealVar frac1("frac1","frac1",7.1345e-01,0.6,0.72);
    RooRealVar frac2("frac2","frac2",1.9309e-01,0.191,0.194);

    RooAddPdf* signal;
    
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
    //RooDataSet combData("combData","combined data",InvariantMass,Index(sample),Import("ALL", *Data_ALL),Import("PASSING", *Data_PASSING));
    
    // C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
    // -----------------------------------------------------------------------------------
    RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
    
    simPdf.addPdf(*model,"ALL");
    simPdf.addPdf(*model_pass,"PASSING");
    
    RooFitResult* fitres = new RooFitResult;
    fitres = simPdf.fitTo(combData, RooFit::Save());
    
    // OUTPUT ARRAY
    RooRealVar** fit_return = new RooRealVar*[2];
    double* output = new double[4];
    
    RooRealVar* yield_ALL = (RooRealVar*) fitres->floatParsFinal().find("n_signal_total");
    RooRealVar* yield_PASS = (RooRealVar*) fitres->floatParsFinal().find("n_signal_total_pass");
    
    output[0] = yield_ALL->getVal();
    output[1] = yield_PASS->getVal();
    
    output[2] = yield_ALL->getError();
    output[3] = yield_PASS->getError();

    
    cout << "VAL = " << output[0] << " +- " << output[2] << endl;
    
    // Some Random Stuff like... whaaaa....
    frame->SetTitle("ALL");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    Data_ALL->plotOn(frame);
    
    //model->paramOn(frame,Layout(0.60,0.90,0.75));
    model->plotOn(frame);
    model->plotOn(frame,RooFit::Components("signal1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model->plotOn(frame,RooFit::Components("signal2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model->plotOn(frame,RooFit::Components("signal3"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
    model->plotOn(frame,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    
    c_all->cd();
    frame->Draw("");
    
    //TESTING THE NEW DATASET
    RooPlot *frame_pass = InvariantMass.frame(RooFit::Title("Invariant Mass"));
    
    c_pass->cd();
    
    frame_pass->SetTitle("PASSING");
    frame_pass->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
    Data_PASSING->plotOn(frame_pass);
    
    //model_pass->paramOn(frame_pass,Layout(0.60,0.70,0.75));
    model_pass->plotOn(frame_pass);
    model_pass->plotOn(frame_pass,RooFit::Components("signal1"),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
    model_pass->plotOn(frame_pass,RooFit::Components("signal2"),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta - 5));
    model_pass->plotOn(frame_pass,RooFit::Components("signal3"),RooFit::LineStyle(kDashed),RooFit::LineColor(kOrange));
    model_pass->plotOn(frame_pass,RooFit::Components("background"),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    
    frame_pass->Draw();
    
    string file     = "Result/";
    string all_pdf  = "_ALL.pdf";
    string pass_pdf = "_PASS.pdf";
    
    if(save)
    {
        c_pass->SaveAs(strcat(strcat(file,condition), all_pdf));
        c_all->SaveAs(strcat(strcat(file, condition),pass_pdf));
    }
        
    // DELETING ALLOCATED MEMORY
    delete file0;
    //delete DataTree; - Deleting TTree deletes TFile
    
    delete Data_ALL;
    delete Data_PASSING;
    //
    delete dh_ALL;
    delete dh_PASSING;
    //
    delete cutvar;
    delete redeuce;
    //
    delete signal;
    //
    delete c_all;
    delete c_pass;
    
    //not sure how to delete these pointers...
    //delete frame;
    //delete frame_pass;
    
    delete model;
    delete model_pass;
    delete fitres;
    
    return output;
}

double* vary_bin(double* values, double factor, int bin_n) // factor € [0, 1]
{
    double *re_distribute = new double[bin_n];
    double *output        = new double[bin_n];
    
    for(int i = 0; i < bin_n/2; i++)
    {
        re_distribute[i]  = ((double)bin_n - (double)i)/(double)bin_n * factor * values[i];
        output[i] = values[i] - re_distribute[i];
        output[bin_n-i-1] = (values[bin_n-i-1] + re_distribute[i]);
    }
    delete[] re_distribute;
    if (bin_n % 2 != 0)
    {
        output[bin_n/2] = (double)(output[bin_n/2 - 1] + output[bin_n/2 + 1])/2.0;
    }
    return output;
}

double* binning(int bin_n, double up_lim)
{
    double* quant_values = new double[bin_n+1];
    
    for (int i = 0; i <= bin_n; i++)
    {
        quant_values[i] = (double)(up_lim - 2.990)/(double)(bin_n);
    }
    quant_values = vary_bin(quant_values, 0.75, bin_n + 1);
    
    quant_values[0] = 2.990;
    for (int i = 1; i < bin_n; i++)
    {
        quant_values[i] += quant_values[i-1];
    }
    quant_values[bin_n] = up_lim;
    return quant_values;
}

string* get_conditions(int bin_n, double up_lim)
{
    //up_lim -= 3.5;
    //int bin_n = 160;
    double* quant_values = binning(bin_n, up_lim);
    
    string* conditions = new string[bin_n];
    
    string quantity = "ProbeMuon_Pt";
    
    for (int i = 0; i < bin_n; i++)
    {
        conditions[i] = quantity + ">" + to_string(quant_values[i]) + " && " + quantity + "<" + to_string(quant_values[i+1]);
    }
    
    delete[] quant_values;
    return conditions;
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

/* // GENERATES EQUALLY SPACED CONDITIONS
 string* get_conditions(int bin_n, double up_lim = 60, string quantity = "ProbeMuon_Pt")
 {
     up_lim -= 3.5;
     //int bin_n = 160;
     double* quant_values = new double[bin_n+1];
     
     string* conditions = new string[bin_n];
     
     for (int i = 0; i <= bin_n; i++)
     {
         quant_values[i] = up_lim*(double)i/(double)bin_n + (double)3.5;
     }
     
     for (int i = 0; i < bin_n; i++)
     {
         conditions[i] = quantity + ">" + to_string(quant_values[i]) + " && " + quantity + "<" + to_string(quant_values[i+1]);
     }
     
     delete[] quant_values;
     return conditions;
 }
*/

TH1F* make_hist(string name, double** values, int qnt, int bin_n, double up_lim, bool DRAW = FALSE)
{
    
    TH1F* hist = new TH1F(name.c_str(), name.c_str(), bin_n, 1, up_lim);
    for (int i = 0; i < bin_n; i++)
    {
        hist->AddBinContent(i, values[i][qnt]);
    }
    if (DRAW)
    {
        TCanvas* xperiment = new TCanvas;
        xperiment->cd();
        hist->Draw();
    }
    return hist;
}

//MANUALLY DEFINE THE BINNING
//SET BIN ERROR ON HISTOGRAM
TH1F* make_hist(string name, double** values, int qnt, int bin_n, Double_t* binning, bool DRAW = FALSE)
{
    //AddBinContent
    //HISTOGRAM NEEDS TO HAVE VARIABLE BINS
   
    TH1F* hist = new TH1F(name.c_str(), name.c_str(), bin_n, binning);

    for (int i = 0; i < bin_n; i++)
    {
        hist->SetBinContent(i, values[i][qnt]);
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

TEfficiency* get_efficiency(TH1F* ALL, TH1F* PASS, TH1F* ALL_ERR = NULL, TH1F* PASS_ERR = NULL)
{
    TFile* pFile = new TFile("Efficiency.root","recreate");
    TEfficiency* pEff = new TEfficiency();
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
    /*
    TFile* pFile = new TFile("Efficiency.root","recreate");
    TEfficiency* pEff = new TEfficiency(*ALL, *PASS);
    pEff->SetDirectory(gDirectory);
    pFile->Write();
    return pEff;
     */
}

void Efficiency()
{
    double bins[] = {2, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.2, 6.4, 6.6, 6.8, 7.3, 7.6, 8.0, 8.5, 9.0, 10.0, 11.0, 13.0, 17.0, 50.0};
    int bin_n = 43;
    
    string* conditions = get_conditions(bin_n, bins);
    double ** yields_n_errs = new double*[bin_n]; // [yield_all, yield_pass, err_all, err_pass]
    
    for (int i = 0; i < bin_n; i++)
    {
        yields_n_errs[i] = doFit(conditions[i], TRUE);
    }

    TH1F *yield_ALL  = make_hist("ALL", yields_n_errs, 0, bin_n, bins);
    TH1F *yield_PASS = make_hist("PASS", yields_n_errs, 1,bin_n, bins);
    
    // saves histograms to .root
    TFile* generatedFile = TFile::Open("Histograms.root","RECREATE");
    yield_ALL->SetDirectory(gDirectory);
    yield_PASS->SetDirectory(gDirectory);
    generatedFile->Write();
    
    TH1F *err_ALL    = new TH1F("err_ALL", "err_ALL", 100, 0, 0.4);
    TH1F *err_Pass   = new TH1F("err_Pass", "err_Pass", 100, 0, 0.4);
    
    get_efficiency(yield_ALL, yield_PASS);
    
    //delete bins;
    delete[] yields_n_errs;
    delete[] conditions;
}
