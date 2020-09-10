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

//    GENERATE T&P_UPSILON_DATA.root FROM NTUPPLE
     
    // This is the story of 2 TTrees and how their love manifested itself into physical existance through their offspring... It all started in a TFile...
void get_root()
{
    TFile *file0  = TFile::Open("DATA/JPsi/Run2011AMuOnia_mergeNtuple.root");

    TTree *DataTree = (TTree*)file0->Get(("tagandprobe/AnalysisTree"));
    TTree *PlotControl = (TTree*)file0->Get(("tagandprobe/PlotControl"));
    
    int N_ENTRIES = DataTree->GetEntries(); // -------------------------- LIMIT AMOUNT OF DATA HERE

    TFile Output("TP_JPSI_DATA.root", "recreate");
    
    TTree UPSILON_DATA("JPSI_DATA", "JPSI_DATA");
    
    Int_t PassingProbeStandAloneMuon;
    Double_t InvariantMass, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi;
    
    DataTree->SetBranchAddress("PassingProbeStandAloneMuon", &PassingProbeStandAloneMuon);
    DataTree->SetBranchAddress("InvariantMass", &InvariantMass);
    
    PlotControl->SetBranchAddress("ProbeMuon_Pt" , &ProbeMuon_Pt);
    PlotControl->SetBranchAddress("ProbeMuon_Eta", &ProbeMuon_Eta);
    PlotControl->SetBranchAddress("ProbeMuon_Phi", &ProbeMuon_Phi);
    
    UPSILON_DATA.Branch("InvariantMass",&InvariantMass);
    UPSILON_DATA.Branch("PassingProbeStandAloneMuon", &PassingProbeStandAloneMuon);
    UPSILON_DATA.Branch("ProbeMuon_Pt" ,&ProbeMuon_Pt);
    UPSILON_DATA.Branch("ProbeMuon_Eta",&ProbeMuon_Eta);
    UPSILON_DATA.Branch("ProbeMuon_Phi",&ProbeMuon_Phi);
    
    // show progression
    string progressFormat = "progress at: %f %"+to_string(strlen(to_string(N_ENTRIES).data()))+"\r";
    
    // fill the tree
    for (Int_t i=0; i<N_ENTRIES; i++)
    {
        printf((progressFormat).data(), ((float)i/(float)N_ENTRIES)*100);
        
        DataTree->GetEntry(i);
        PlotControl->GetEntry(i);
        UPSILON_DATA.Fill();
    }
    UPSILON_DATA.Write();
}
