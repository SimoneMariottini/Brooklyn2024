#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>
#include "AnaTools.h"


#define NCHANNELS 16

using namespace std;

void PlotSignificance(){

    TString inname("run_190424_trgch0-15_thr15_gate80.root");
    TString path("/Users/simonemariottini/Documents/Università/Magistrale/Lab2/Measurements/runs/");


    TFile* f = new TFile(path + inname, "READ");

    f->cd();

    TDirectory* dir;
    dir = f->GetDirectory("Hist_Charge");

    if(dir){
        dir->cd();

        for(int i = 0; i < NCHANNELS; i++){
            auto c1 = new TCanvas(Form("Channel_%i", i), Form("Channel_%i", i), 1200, 800);
            c1->cd();
            TH1D * h1 = (TH1D*)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));

            h1->Draw();
        
            TF1* langausgaus = h1->GetFunction("langausgaus");
            TF1* pedgaus = new TF1(Form("pedgaus_%i", i), "[0]*TMath::Gaus(x,[1],[2])", -0.01, 0.05);
            TF1* langaus = new TF1(Form("langaus_%i", i), AnaTools::lanGausFun, -0.01, 0.05, 4);

            double* param = langausgaus->GetParameters();
            langaus->SetParameters(param);
            pedgaus->SetParameters(param + 4);

            TH1D* sig = new TH1D(Form("Significance_channel_%i", i), Form("Significance channel %i", i), 500, -0.01, 0.05);

            for(int j = 1; j <= 500; j++){
                sig->SetBinContent(j, langaus->Integral(sig->GetBinCenter(j), 0.05, 1.0E-1)/TMath::Sqrt(langausgaus->Integral(sig->GetBinCenter(j), 0.05, 1.0E-1)));
            }

            double max = sig->GetBinCenter(sig->GetMaximumBin());

            //your code 
            TLine *v_line= new TLine(max,0,max,200); //declare the vertical line 
            //Set line attributes 
            v_line->SetLineColor(kBlue);
            v_line->SetLineWidth(2);
            v_line->SetLineStyle(kDashed);
            //draw the line
            v_line->Draw("same");
        }
    }
    else{
        cout << "Couldn't find directory" << endl;
    }
}