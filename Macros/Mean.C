#include <TNtupleD.h>

using namespace std;

#define NLAYERS 16

void Mean(){

    TFile* f = new TFile("./Runs/geant_energy.root", "READ");
    f->cd();

    TNtupleD* ntuple = (TNtupleD*)gDirectory->Get("Brooklyn2024");

    ROOT::RDataFrame* dataframe = new ROOT::RDataFrame(*ntuple);

    for (int i = 0; i < NLAYERS; i++){
        cout << "Channel " << i << ":" << endl; 
        cout << "E Mean = " << dataframe->Mean(Form("E_%i", i)).GetValue() << endl;
        cout << "L Mean = " << dataframe->Mean(Form("L_%i", i)).GetValue() << endl;
    }
}