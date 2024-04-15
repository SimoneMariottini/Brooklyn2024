void prova() {

    TFile *file = new TFile("hist.root", "READ");
    file->cd();
    TH1D *h = (TH1D *) gDirectory->Get("histo");
    h->Draw();

    TF1 *fit = new  TF1("fit", "gaus(0)", -5., 5.);

    fit->SetParameter(0, h->GetBinContent(h->GetMaximumBin()));
    fit->SetParameter(1, h->GetBinCenter(h->GetMaximumBin()));
    fit->SetParameter(2, h->GetStdDev());
    fit->SetParLimits(1, -2., 2.);
    h->Fit(fit, "R");

}