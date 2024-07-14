#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include <TString.h>
#include <TFile.h>
#include "../Header/AnaTools.h"
#include <iostream>

#define NCHANNELS 16

using namespace std;

void NewFitChargeLandau(){

    TString inname("run_190424_trgch0-15_thr15_gate80.root");
    TString path("./Runs/");

    TFile* f = new TFile(path + inname, "READ");

    f->cd();

    AnaTools MyAnaTools(f); 

    TF1 *langausgaus = new TF1("langausgaus", AnaTools::lanGausPlusGausFun, -0.01, 0.05, 7);

    double par[7] = {0.0006, 0.001, 0.54, 0.015, 0.65, -0.0008, 0.0004};

    langausgaus->SetNpx(10000);
    langausgaus->SetParameters(par);
    langausgaus->SetParLimits(3, 0.01, 10);
    langausgaus->SetParNames("Width", "MP", "Area", "G Sigma", "Ped Constant", "Ped Mean", "Ped Sigma");
    
    TCanvas* c1 = new TCanvas();
    c1->cd();

    TH1D* histo = new TH1D(*(MyAnaTools.GetChargeHistogram(5)));

    auto FitResult = histo->Fit(langausgaus, "", "", -0.01, 0.5);
    
    histo->Draw();
}
