
```cpp
double* doFit(string condition, string MuonID_str, bool save = TRUE) // RETURNS ARRAY WITH [yield_all, yield_pass, err_all, err_pass]    ->   OUTPUT ARRAY
```

void compare_efficiency()

void get_TTree_from_ntupple()

TEfficiency* get_efficiency(TH1F* ALL, TH1F* PASS)

string* get_conditions(int bin_n, double* bins, string quantity)

TH1F* make_hist(string name, double** values, int qnt, int bin_n, Double_t* binning, bool IsDataMc, bool DRAW = FALSE)

double* McYield(string condition)

void change_bin(string condition, string hist_file = "Histograms_Run2011.root")
