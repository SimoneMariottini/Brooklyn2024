#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "../Header/AnaTools.h"
#include <TFile.h>

#define NCHANNELS 16

void EvaluateTime(){

TFile *f = new TFile("./Runs/run_190424_trgch0-15_thr15_gate80.root", "READ");

AnaTools MyAnatools = AnaTools(f, "info.root"); //legge il file 


}