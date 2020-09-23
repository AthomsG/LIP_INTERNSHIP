
```cpp
double* doFit(string condition, string MuonID_str, bool save = TRUE) // RETURNS ARRAY WITH [yield_all, yield_pass, err_all, err_pass]    ->   OUTPUT ARRAY
```

```cpp
void compare_efficiency()
```

```cpp
void get_TTree_from_ntupple()
```

```cpp
TEfficiency* get_efficiency(TH1F* ALL, TH1F* PASS)
```

```cpp
string* get_conditions(int bin_n, double* bins, string quantity)
```

```cpp
TH1F* make_hist(string name, double** values, int qnt, int bin_n, Double_t* binning, bool IsDataMc, bool DRAW = FALSE)
```

```cpp
double* McYield(string condition)
```

```cpp
void change_bin(string condition, string hist_file = "Histograms_Run2011.root")
```
