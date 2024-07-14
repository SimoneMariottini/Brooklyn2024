#ifdef __CLING__
#pragma cling optimize(0)
#endif

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

double langausfunction(double* x, double* par);

void DrawChargeFitResult(){

    //TString inname("run_190424_trgch0-15_thr15_gate80.root");
    TString path("./FitResults/");

        for(int i = 1; i < NCHANNELS - 1; i++){
            TFile* f = new TFile(path + Form("fit_landau_%i.root", i), "READ");

            f->cd();
            auto c1 = new TCanvas();
            c1->cd();
            TH1D * h1 = (TH1D*)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));

            h1->Draw();
        
            TF1* langausgaus = h1->GetFunction("langausgaus");
            TF1* pedgaus = new TF1(Form("pedgaus_%i", i), "[0]*TMath::Gaus(x,[1],[2], true)", -0.01, 0.05);
            TF1* langaus = new TF1(Form("langaus_%i", i), langausfunction, -0.01, 0.05, 6);

            pedgaus->SetLineColor(2);
            langaus->SetLineColor(3);
            langausgaus->SetLineColor(4);
            h1->SetLineColor(1);

            double* param = langausgaus->GetParameters();

            double par[6] = {param[0], param[1] - param[5], param[2], param[3], param[6], param[5]};

            langaus->SetParameters(par);
            pedgaus->SetParameters(param + 4);

            langaus->SetNpx(1000);
            pedgaus->SetNpx(1000);

            langausgaus->SetLineWidth(2);
            langaus->SetLineWidth(2);
            pedgaus->SetLineWidth(2);

            langausgaus->SetParNames("Width", "Most Prob. Value", "Area", "Resolution Factor", "Ped. Area", "Ped. #mu", "Ped. #sigma");

            pedgaus->Draw("same");
            langaus->Draw("same");

            TLegend* legend = new TLegend(0.3, 0.2, "");
            legend->AddEntry(Form("langaus_%i", i), "Convoluted Landau");
            legend->AddEntry(Form("pedgaus_%i", i), "Pedestal Gaussian");
            legend->AddEntry(langausgaus, "Total Distribution");
            legend->AddEntry(h1, "Data");
            legend->Draw();

            gStyle->SetOptFit(11111111);
        }
    }

double langausfunction(double* x, double* par){
    double x_shift = x[0] - par[5];
    return AnaTools::lanGausFun(&x_shift , par);
}