#include "../Header/AnaTools.h"

void FitChargePoisson(){

    int npar = 6;
    TString infile = "";
    TFile *f = new TFile(infile, "READ");
    f->cd();

    AnaTools *data = new AnaTools(f);

    TCanvas *c1 = new TCanvas();

    TF1 * fun = new TF1("PoisGaus", AnaTools::poisGausFun, -0.01, 0.05, npar);


}