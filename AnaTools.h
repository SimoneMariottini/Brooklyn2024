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
  AnaTools(TFile *f, Event *myEvent, double cf, double th);
  virtual ~AnaTools();

  //Data Analysis Methods
  void BookingHistograms();

  void BookWaveform();
  void BookPersistence();
  void BookCharge();
  void BookToF();

  void Process(int);
  
  void Clear();
  double ComputeTimeCFD(Waveform*, double);
  double ComputeTimeFT(Waveform*, double);
  double ComputeCharge(Waveform*);

  void LoadPedestal(string inname);

 private:
  
  double cf_;
  double th_;

  TFile *outfile_;
  Event *event_;
  int nev_;


  TH1D* h_c_vector_[NCHANNELS];
  TH1D *h_c_tot_;

  TH1D *h_TOF_cfm_[NCHANNELS];
  TH1D *h_TOF_ft_[NCHANNELS];

  TH1D *h_wave_vector_[20][NCHANNELS];

  TH2D *persistence_vector_[NCHANNELS];
  
  
  // double ped_mean[16];
  // double ped_sigma[16];



};

#endif
