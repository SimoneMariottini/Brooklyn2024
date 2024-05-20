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
    TString path("./");

    TFile* f = new TFile(path+inname, "READ");

    AnaTools myAnaTools = AnaTools(f); 

    double* maxSig;
    maxSig = myAnaTools.EvaluateMaxSignificanceBinCenter();

    for(int i = 0; i< NCHANNELS; i++){
        myAnaTools.SetCutoff(i, maxSig[i]);
            /*
            //your code 
            TLine *v_line= new TLine(maxSig[i],0,maxSig[i],200); //declare the vertical line 
            //Set line attributes 
            v_line->SetLineColorkBlue);
            v_line->SetLineWidth(2);
            v_line->SetLineStyle(kDashed);
            //draw the line
            v_line->Draw("same");
            */
    }
    myAnaTools.SaveInfo("test_info.root", "RECREATE");
}