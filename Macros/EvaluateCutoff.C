#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include "../Header/AnaTools.h"

#define NCHANNELS 16

using namespace std;
 
void EvaluateCutoff(){

    //TString inname("run_190424_trgch0-15_thr15_gate80.root");
    TString path("./FitResults/");

    double cutoff[NCHANNELS] = {0.};

    AnaTools* myAna = new AnaTools("info.root");

    for(int i = 1; i < NCHANNELS - 1; i++){
        TFile* f = new TFile(path + Form("fit_landau_%i.root", i), "READ");

        f->cd();

        TH1D * h1 = (TH1D*)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));
        
        TF1* langausgaus = h1->GetFunction("langausgaus");

        cutoff[i] = langausgaus->GetMinimumX(langausgaus->GetParameter(5), langausgaus->GetParameter(1));

        myAna->SetCutoff(i, cutoff[i]);
    }
 
    myAna->SaveInfo("cutoff");
}
