#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include "AnaTools.h"

#define NCHANNELS 16

using namespace std;

void DrawChargeFitResult(){

    TString inname("run_190424_trgch0-15_thr15_gate80.root");
    TString path("/Users/simonemariottini/Documents/Università/Magistrale/Lab2/Measurements/runs/");


    TFile* f = new TFile(path + inname, "READ");

    f->cd();

    TDirectory* dir;
    dir = f->GetDirectory("Hist_Charge");

    if(dir){
        dir->cd();
    
        gStyle->SetOptFit(1111111);

        for(int i = 0; i < NCHANNELS; i++){
           auto c1 = new TCanvas(Form("Channel_%i", i), Form("Channel_%i", i), 1200, 800);
           c1->cd();
           TH1D * h1 = (TH1D*)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));

          h1->Draw();
        
          TF1* langausgaus = h1->GetFunction("langausgaus");
          TF1* pedgaus = new TF1(Form("pedgaus_%i", i), "[0]*TMath::Gaus(x,[1],[2])", -0.01, 0.05);
          TF1* langaus = new TF1(Form("langaus_%i", i), AnaTools::lanGausFun, -0.01, 0.05, 4);

          pedgaus->SetLineColor(kBlue);
            langaus->SetLineColor(kGreen + 3);

            double* param = langausgaus->GetParameters();
            langaus->SetParameters(param);
            pedgaus->SetParameters(param + 4);

            langaus->SetNpx(1000);
            pedgaus->SetNpx(1000);

            langaus->SetLineWidth(4);
            pedgaus->SetLineWidth(4);

            pedgaus->Draw("same");
            langaus->Draw("same");
        }
    }
}