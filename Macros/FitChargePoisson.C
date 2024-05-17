#include "../Header/AnaTools.h"

void FitChargePoisson(){

    /*TString infile = "";
    TFile *f = new TFile(infile, "READ");
    f->cd();

    AnaTools *data = new AnaTools(f);*/

    TCanvas *c1 = new TCanvas();

    TF1 * fun = new TF1("PoisGaus", AnaTools::poisGausFun, -0.01, 0.05, 6);

    double par[6] = {-0.0018, 0.0005, 0.0018, 0.0008, 2, 1};

    fun->SetParameters(par);

    fun->SetNpx(10000);

    fun->Draw();


}