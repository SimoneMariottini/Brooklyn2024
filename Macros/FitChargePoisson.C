#include "../Header/AnaTools.h"
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>

#define NCHANNELS 16

using namespace std;

void FitChargePoisson(){

    double invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)

    TString path("/Users/robertomattacchione/Desktop/lab2/");
    TString inname("run_160524_trgledonch15_lowlight_long.root");
    TFile *f = new TFile(path + inname, "READ");
    f->cd();

    AnaTools *data = new AnaTools(f);

    TDirectory* dir;
    dir = f->GetDirectory("Hist_Charge");

    TF1 * fun = new TF1("PoisGaus", AnaTools::poisGausFun, -0.01, 0.05, 6);

    //double par[6] = {-0.0018, 0.0005, 0.0018, 0.0008, 2, 1};
    double par[6] = {0, 0.0005, 0.0017, 0.0005, 0, 12};

    fun->SetParameters(par);
    //fun->SetParLimits(2,0.0017,0.0018);
    fun->SetNpx(10000);
    fun->SetParNames("Pedestal #mu", "Pedestal #sigma", "Photon charge #mu", "Photon charge #sigma", "#lambda", "Integral");

    if(dir){
        dir->cd();
    
        gStyle->SetOptFit(11111111);

        for (int i = 0; i < NCHANNELS; i++)
        {   
            //int i = 12;

            auto c1 = new TCanvas(Form("Channel_%i", i), Form("Channel_%i", i), 1200, 800);
    
            TH1D* h1 = (TH1D*)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));

            h1->Fit("PoisGaus","","same",-0.01,0.05);

            h1->Draw();
            
            TF1* func = h1->GetFunction("PoisGaus");

            TF1* g0 = new TF1(Form("gaus0_%i", i),"TMath::Exp(-[2])*[3]*TMath::Gaus(x,[0],[1],1)",-0.01,0.05);

            //sigma1 = TMath::Sqrt(TMath::Power([1], 2) + TMath::Power([3], 2));
            TF1* g1 = new TF1(Form("gaus1_%i", i),"[4]*TMath::Exp(-[4])*[5]*TMath::Gaus(x,[0]+[2],TMath::Sqrt(TMath::Power([1], 2) + TMath::Power([3], 2)),1)",-0.01,0.05);

            //sigma2 = TMath::Sqrt(TMath::Power([1], 2) + TMath::Power([3], 2) * 2);
            TF1* g2 = new TF1(Form("gaus2_%i", i),"[4]*[4]*(1/2.)*TMath::Exp(-[4])*[5]*TMath::Gaus(x,[0]+2*[2],TMath::Sqrt(TMath::Power([1], 2) + TMath::Power([3], 2) * 2),1)",-0.01,0.05);

            //sigma2 = TMath::Sqrt(TMath::Power([1], 2) + TMath::Power([3], 2) * 2);
            TF1* g3 = new TF1(Form("gaus3_%i", i),"[4]*[4]*[4]*(1/6.)*TMath::Exp(-[4])*[5]*TMath::Gaus(x,[0]+3*[2],TMath::Sqrt(TMath::Power([1], 2) + TMath::Power([3], 2) * 3),1)",-0.01,0.05);

            double* param = func->GetParameters();

            g0->SetParameter(0, param[0]);
            g0->SetParameter(1, param[1]);
            g0->SetParameter(2, param[4]);
            g0->SetParameter(3, param[5]);
            g1->SetParameters(param);
            g2->SetParameters(param);
            g3->SetParameters(param);

            g0->SetLineColorAlpha(4,0.8);
            g1->SetLineColorAlpha(2,0.8);
            g2->SetLineColorAlpha(3,0.8);
            g3->SetLineColorAlpha(6,0.8);

            g0->SetNpx(10000);
            g1->SetNpx(10000);
            g2->SetNpx(10000);
            g3->SetNpx(10000);
            

            g0->SetLineWidth(2);
            g1->SetLineWidth(2);
            g2->SetLineWidth(2);
            g3->SetLineWidth(2);

            g0->Draw("same");
            g1->Draw("same");
            g2->Draw("same");
            g3->Draw("same");

            fun->SetParameters(par);
        }
    }   
}
