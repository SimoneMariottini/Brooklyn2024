#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include <TString.h>
#include <TFile.h>
#include "AnaTools.h"
#include <iostream>

#define NCHANNELS 16

using namespace std;

void FitChargeMacro(){

    TString inname("run_190424_trgch0-15_thr15_gate80.root");
    TString path("/Users/simonemariottini/Documents/Università/Magistrale/Lab2/Measurements/runs/");

    TFile* f = new TFile(path + inname, "UPDATE");

    f->cd();

    AnaTools MyAnaTools(f); 

    double* efficiency = MyAnaTools.EvaluateEfficiency();

    for(int i = 0; i< NCHANNELS; i++){
    cout << "Channel " << i << " efficiency = " << efficiency[i] << "+/-" << efficiency[i + NCHANNELS] << endl;
    }

    f->Write();
    f->Close();
}
