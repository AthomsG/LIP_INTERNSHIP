# Calculating Efficiencies using Tag & Probe

> Tag &amp; probe efficiency fitting method

## Setup

This project was developed using [ROOT](https://root.cern.ch/root/html534/guides/users-guide/InstallandBuild.html), made available by CERN and the following Datasets:
* [1] [Run2011AMuOnia_mergeNtuple.root](https://drive.google.com/drive/u/0/folders/1Nu9Al7SV1F60TMFxKZVBIMvgEWAdzida)
* [2] [JPsiToMuMu_mergeMCNtuple.root](https://drive.google.com/drive/u/0/folders/1Nu9Al7SV1F60TMFxKZVBIMvgEWAdzida)

From these two datasets, a `.root` file was generated for each MuonId (i.e *Standalone*, *Tracker* or *Global*) using the `get_root.C` and then stored on the `\Data` folder, following the respective hierarchy.

## Fitting Method




## WorkFlow

`Efficiency.C` is given as an example of how to use the fitting method to calculate an efficiency.
First the user has to manually define the bins in which the quantity being studied (i.e. pT, Eta, Phi) will be divided

## Preferences

You can change the method to estimate signal region by modifying `Muon.setMethod(1)` line by choosing 1 (estimate by FWHM of histograms) or 2 (estimate by FWHM of fitting):

```cpp
Muon.setMethod(1);
```

Change this line to specify the ntupple you are analysing by choosing 0 (old ntupple), 1 (run 2011 ntupple) or 2 (monte carlo ntupple):

```cpp
int useNewData = 1;
```

## Running

Go on your folder where the file code is downloaded and run:

```sh
$ cd main
$ root -l -n
root[0] .L macro.cpp+
root[1] macro()
```

## Output
Output images are stored in the `result` folder.
