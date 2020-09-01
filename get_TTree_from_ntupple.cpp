
//    GENERATE T&P_UPSILON_DATA.root FROM NTUPPLE
     
    // This is the story of 2 TTrees and how their love manifested itself into physical existance through their offspring... It all started in a TFile...
{
    TFile *file0  = TFile::Open("Run2011_upsilon.root");

    TTree *DataTree = (TTree*)file0->Get(("tagandprobe/AnalysisTree"));
    TTree *PlotControl = (TTree*)file0->Get(("tagandprobe/PlotControl"));
    
    int N_ENTRIES = DataTree->GetEntries(); // -------------------------- LIMIT AMOUNT OF DATA HERE

    TFile Output("T&P_UPSILON_DATA.root", "recreate");
    
    TTree UPSILON_DATA("UPSILON_DATA", "UPSILON_DATA");
    
    Int_t PassingProbeTrackingMuon;
    Double_t InvariantMass, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi;
    
    DataTree->SetBranchAddress("PassingProbeTrackingMuon", &PassingProbeTrackingMuon);
    DataTree->SetBranchAddress("InvariantMass", &InvariantMass);
    
    PlotControl->SetBranchAddress("ProbeMuon_Pt" , &ProbeMuon_Pt);
    PlotControl->SetBranchAddress("ProbeMuon_Eta", &ProbeMuon_Eta);
    PlotControl->SetBranchAddress("ProbeMuon_Phi", &ProbeMuon_Phi);
    
    UPSILON_DATA.Branch("InvariantMass",&InvariantMass);
    UPSILON_DATA.Branch("PassingProbeTrackingMuon", &PassingProbeTrackingMuon);
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