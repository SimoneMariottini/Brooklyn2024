#include <TTree.h>

void DrawFromTree(){
    

    TFile* f = new TFile("./Runs/geant_energy.root", "READ");
    f->cd();
    TH1D* hist = new TH1D("histo", "Simulated Energy Deposit, Channel 5", 500, 0., 0.35);

    hist->SetXTitle("#DeltaE [MeV]");
    hist->SetYTitle("Counts [#]");

    TTree* tree = (TTree*)gDirectory->Get("Brooklyn2024");

    tree->Draw("E_5>>histo");

    hist->Draw();
}