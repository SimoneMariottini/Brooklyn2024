#include "../Header/AnaTools.h"

void MergeResult(){
    TFile* final = new TFile("Landau_fit_results.root", "UPDATE");
    
    TFile* result = new TFile("./FitResults/fit_landau_1.root", "READ");

    final->cd();

    TH1D* hist = new TH1D(*(TH1D*)(result->Get("Hist_Charge_Channel_1")));

    TF1 *func = (TF1*)hist->GetFunction("langausgaus");
    func->Print("");
    //hist->Fit(func, "", "same", -0.01, 0.05);

    hist->Draw();




}