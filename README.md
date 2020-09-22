# Calculating Efficiencies using Tag & Probe

> Tag &amp; probe efficiency fitting method

## Setup

This project was developed using [ROOT](https://root.cern.ch/root/html534/guides/users-guide/InstallandBuild.html), made available by CERN and the following Datasets:
* [1] [Run2011AMuOnia_mergeNtuple.root](https://drive.google.com/drive/u/0/folders/1Nu9Al7SV1F60TMFxKZVBIMvgEWAdzida)
* [2] [JPsiToMuMu_mergeMCNtuple.root](https://drive.google.com/drive/u/0/folders/1Nu9Al7SV1F60TMFxKZVBIMvgEWAdzida)

From these two datasets, a `.root` file was generated for each MuonId (i.e *Standalone*, *Tracker* or *Global*) using the `get_root.C` and then stored on the `\Data` folder, following the respective hierarchy.

## Fitting Method

![Alt text](images/esquema.png){ width: 200px; height:100px}


## WorkFlow

`Efficiency.C` is given as an example of how to use the fitting method to calculate an efficiency. It follows as such:
1. The user has to manually define the bins in which the quantity being studied (i.e. pT, Eta, Phi) will be divided;
2. Generate conditions (that divide the dataset into the defined binned intervals) using ```get_conditions```;
3. Create a loop that fits the invariant mass for each bin using ```doFit```when the dataset consists of real data and ```McYield``` for the Monte Carlo dataset.

## Running

On this repository, do:

```sh
root -l -b -q Efficiency.C
```

## Output
Output images are stored in the `/result` folder.
