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

  AnaTools(TFile *f, Event *myEvent, double cf, double th); //Normal constructor

  AnaTools(TFile *f); //"Copy from .root file" constructor

  AnaTools(TFile *f, Event *myEvent, double cf, double th, TString infoFile); //Add path to detector_info.root file to save/write to
  AnaTools(TFile *f, TString infoFile); //Add path to detector_info.root file to save/write to

  virtual ~AnaTools();

  //Data Analysis Methods
  void BookWaveform();
  void BookPersistence();
  void BookCharge();
  void BookToF();

  void Process(int);
  
  void Clear();
  double ComputeTimeCFD(Waveform*, double);
  double ComputeTimeFT(Waveform*, double);
  double ComputeCharge(Waveform*);

  TFitResultPtr* FitCharge();

  static double lanGausFun(double *x, double *par);

  static double lanGausPlusGausFun(double *x, double *par);

  double* EvaluateEfficiency();

  void LoadInfo(TString);

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
  
  bool bookings_[4] = {0}; //0: Waveform, 1: Persistence, 2: Charge, 3: ToF

  //Info on system characterization

  TString infoFile_ = "";

  double cutoff_[NCHANNELS];
};

#endif
