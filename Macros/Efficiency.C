#include "../Header/AnaTools.h"
#include <TMatrix.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

double langausfunction(double* x, double* par);

#define NCHANNELS 16

using namespace std;

void Efficiency(){
     TString path("./FitResults/");

        for(int i = 1; i < NCHANNELS - 1; i++){
            TFile* f = new TFile(path + Form("fit_landau_%i.root", i), "READ");

            f->cd();

            TH1D * h1 = (TH1D*)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));
        
            TF1* langausgaus = h1->GetFunction("langausgaus");
            TF1* pedgaus = new TF1(Form("pedgaus_%i", i), "[0]*TMath::Gaus(x,[1],[2], true)", -0.01, 0.05);
            TF1* langaus = new TF1(Form("langaus_%i", i), langausfunction, -0.01, 0.05, 6);

            double* param = langausgaus->GetParameters();

            double par[6] = {param[0], param[1] - param[5], param[2], param[3], param[6], param[5]};

            langaus->SetParameters(par);
            pedgaus->SetParameters(param + 4);

            langaus->SetNpx(1000);
            pedgaus->SetNpx(1000);

            double pedArea = pedgaus->Integral(-0.01, 0.05);
            //double sigmaPed = langausgaus->IntegralError(-0.01, 0.05);
            double langausArea = langaus->Integral(-0.01, 0.05);
            //double sigmaLangaus = langaus->IntegralError(-0.01, 0.05);

            double efficiency = langausArea / (langausArea + pedArea);
            //double efficiency_error = 1 / TMath::Power((langausArea + pedArea), 2) * TMath::Sqrt(TMath::Power(pedArea, 2) * TMath::Power(sigmaLangaus, 2) + TMath::Power(langausArea, 2) * TMath::Power(sigmaPed, 2));
            
            cout << "Channel " << i << " Efficiency = " << efficiency << endl;
        }

}

double langausfunction(double* x, double* par){
    double x_shift = x[0] - par[5];
    double new_par[6] = {par[0], par[1] - par[5], par[2], par[3], par[6], par[5]};
    return AnaTools::lanGausFun(&x_shift , par);
}