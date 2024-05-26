#include "../Header/AnaTools.h"

#define NCHANNELS 16

void FitTimeGaussian(){

    TString path("./Runs/");
    TString inname("run_190424_trgch0-15_thr15_gate80.root");
    TFile *f = new TFile(path + inname, "READ");
    f->cd();

    AnaTools *data = new AnaTools(f);

 
    TDirectory* dir;
    dir = f->GetDirectory("./Hist_Time/Hist_Time_Method_1");

    //TF1 * fun = new TF1("Gaus", AnaTools::GausFun, -0.01, 0.05, 6);
    TF1 *fun = new TF1("Gaus", "[2]*TMath::Gaus(x,[0],[1])");

    double par[2] = {-0.8, 0.8};

    fun->SetParameters(par);
    fun->SetNpx(10000);
    fun->SetParNames("#mu", "#sigma");

    if(dir){
        dir->cd();
    
        gStyle->SetOptFit(11111111);
        

        for (int i = 0; i < NCHANNELS; i++)
        {
            auto c1 = new TCanvas(Form("Channel_%i", i), Form("Channel_%i", i), 1200, 800);
    
            TH1D * h1 = (TH1D*)gDirectory->Get(Form("Hist_Time_Channel_%i", i));

            h1->Fit("Gaus","","same",-0.01,0.05);

            h1->Draw();

            fun->SetParameters(par);
        }
    }   
}
