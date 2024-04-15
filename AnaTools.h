#ifndef ANATOOLS_H
#define ANATOOLS_H

#include <stdint.h>
#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include "Waveform.h"
#include "Event.h"
#include <fstream>
#include <sstream>

using namespace std;

class AnaTools{
  
 public:
  
  //Constructor and Destructor
  AnaTools(TFile *f, Event *myEvent, double cf_, double th_);
  virtual ~AnaTools();

  //Data Analysis Methods
  void BookingHistograms();
  void Process(int);
  
  void Clear();
  double ComputeTimeCFD(Waveform*, double);
  double ComputeTimeFT(Waveform*, double);
  double ComputeCharge(Waveform*);

  void LoadPedestal(string inname);

  void Draw();

 private:
  
  double cf;
  double th;
  TFile *outfile;
  Event *event;
  int nev;


  TH1D* hc_vector[NCHANNELS];
  TH1D *hctot;
  TH1D *hTOF_cfm[NCHANNELS];
  TH1D *hTOF_ft[NCHANNELS];
  TH1D *hist_vector[20][NCHANNELS];

  TH2D *persistence_vector[NCHANNELS];
  
  
  // double ped_mean[16];
  // double ped_sigma[16];



};

#endif
