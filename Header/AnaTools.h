#ifndef ANATOOLS_H
#define ANATOOLS_H

#include <stdint.h>
#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include "../Header/Waveform.h"
#include "../Header/Event.h"
#include <fstream>
#include <sstream>

using namespace std;

class AnaTools{
  
 public:
  
  //Constructor and Destructor
  AnaTools();
  AnaTools(TString infoFile){infoFile_ = infoFile;}

  AnaTools(TFile *f, Event *myEvent, double cf, double th); //Normal constructor

  AnaTools(TFile *f); //"Copy from .root file" constructor

  AnaTools(TFile *f, Event *myEvent, double cf, double th, TString infoFile); //Add path to detector_info.root file to save/write to
  AnaTools(TFile *f, TString infoFile); //Add path to detector_info.root file to save/write to

  virtual ~AnaTools();

  //Data Analysis Methods
  void BookWaveform();
  void BookPersistence();

  void BookCharge();
  void BookCharge(const double&, const double&);

  void BookToF();
  void BookTime();

  void Process(int);
  void Clear();
  
  double ComputeTimeCFD(Waveform*, double);
  double ComputeTimeFT(Waveform*, double);
  double ComputeCharge(Waveform*);

  TFitResultPtr* FitCharge();

  static double lanGausFun(double *x, double *par);

  static double lanGausPlusGausFun(double *x, double *par);

  static double poisGausFun(double *x, double *par);

  static double stepFun(double *x);

  static double GausFun(double *x, double *par);

  double* EvaluateEfficiency();

  double* EvaluateMaxSignificanceBinCenter();

  double* EvaluateToF();

  void LoadInfo(TString);
  void LoadInfo(){LoadInfo(infoFile_);};
  void SaveInfo(TString infoType, TString infoFile = "", TString mode = "UPDATE");

  //Getters:
  const TString GetInfoFile() {return infoFile_;};
  const double GetCutoff(int i) {return cutoff_[i];};
  const double GetGain(int i) {return gain_[i][0];};
  const double GetGainError(int i) {return gain_[i][1];};
  const double GetPedestal(int i) {return pedestal_[i];};

  const TH1D* GetChargeHistogram(const int& i){return h_c_vector_[i];}
  const TH1D* GetChargeHistogram(const TString name){ if(name.Contains("tot")) return h_c_tot_; else return nullptr;}

  //Setters:
  void SetCutoff(int i, const double& x);
  void SetEfficiency(int i, const double &x);
  void SetGain(int i, const double &x, const double &sigma_x);
  void SetPedestal(int i, const double &x);

  void SetInfoFile(TString);

  void SetChargeRange(const double&, const double&);

 private:
  
  double cf_;
  double th_;

  TFile *outfile_;
  Event *event_;
  int nev_;

  TH1D* h_c_vector_[NCHANNELS];
  TH1D *h_c_tot_;

  double h_charge_range_[2] = {-0.01, 0.05};

  TH1D *h_TOF_cfm_[NCHANNELS];
  TH1D *h_TOF_ft_[NCHANNELS];

  TH1D *h_time_;
  TH1D *h_time_vector_[2][NCHANNELS];

  TH1D *h_wave_vector_[20][NCHANNELS];

  TH2D *persistence_vector_[NCHANNELS];

  bool bookings_[5] = {0}; //0: Waveform, 1: Persistence, 2: Charge, 3: ToF, 4: Time

  //Info on system characterization
  TString infoFile_ = "info_file.root";
  double cutoff_[NCHANNELS] = {0};
  double efficiency_[NCHANNELS] = {0};
  double gain_[NCHANNELS][2] = {{0}};
  double pedestal_[NCHANNELS] = {0};
};

#endif
