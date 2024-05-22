#include "../Header/AnaTools.h"

#define NCHANNELS 16

void FitChargePoisson(){

    TString path("/Users/robertomattacchione/Desktop/lab2/");
    TString inname("run_170524_trgledonch0_lowlight_long.root");
    TFile *f = new TFile(path + inname, "READ");
    f->cd();

    AnaTools *data = new AnaTools(f);

    double par[6] = {-0.0001, 0.0005, 0.0017, 0.0005, 1, 11.9}; //Function's initial parameters

    TDirectory* dir;
    dir = f->GetDirectory("Hist_Charge");

    TF1 * fun = new TF1("PoisGaus", AnaTools::poisGausFun, -0.01, 0.05, 6);

    double par[6] = {-0.0018, 0.0005, 0.0018, 0.0008, 2, 1};
    //double par[6] = {0.0015, 0.0005, 0.0022, 0.0005, 2, 12};

    fun->SetParameters(par);
    fun->SetParLimits(2,0.0017,0.0018);
    fun->SetNpx(10000);
    fun->SetParNames("Pedestal #mu", "Pedestal #sigma", "Photon charge #mu", "Photon charge #sigma", "#lambda", "Integral");

    if(dir){
        dir->cd();
    
        gStyle->SetOptFit(11111111);

        for (int i = 0; i < NCHANNELS; i++)
        { 
            auto c1 = new TCanvas(Form("Channel_%i", i), Form("Channel_%i", i), 1200, 800);
    
            TH1D * h1 = (TH1D*)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));

            h1->Fit("PoisGaus","","same",-0.01,0.05);

            h1->Draw();

            fun->SetParameters(par);
        }
    }   
}
