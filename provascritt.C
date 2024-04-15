void provascritt() {

    TFile *file = new TFile("hist.root", "RECREATE");
    file->cd();

    TH1D *h = new TH1D("histo", "titolo", 100, -5., 5.);
    h->FillRandom("gaus", 10000);

    file->Write();
    file->Close();
}