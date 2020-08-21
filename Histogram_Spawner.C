#include "RooFitResult.h"

//EFFICIENCY TEST

//PROCURAR TEFFICIENCY
void Histogram_Spawner()
{
    const char *files[3] = {"data_histoall.root",
    "DATA_RAW/Run2011AMuOnia_mergeNtuple.root",
    "JPsiToMuMu_mergeMCNtuple.root"};
    
    //Which file of files (variable above) should use
    int useFile = 1;
    
    //Open and read files
    TFile *file0  = TFile::Open(("../" + string(files[useFile])).data());
    TTree *TreePC = (TTree*)file0->Get(("PlotControl").data());
    TTree *TreeAT = (TTree*)file0->Get(("AnalysisTree").data());
    cout << "Using \"" << files[useFile] << "\" ntupple" << endl;
    
}
