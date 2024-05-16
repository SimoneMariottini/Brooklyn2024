#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <stdint.h>
#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>


#define NSAMPLING 1024 
#define SAMPLINGPERIOD 0.3125
#define NCHANNELS 16
#define ECHARGE 1.6e-19

using namespace std;

class Waveform{
  
 public:
  
  //COSTRUTTORE E DISTRUTTORE
  Waveform(double *vt, double *va);
  virtual ~Waveform();
  
  //GETTERS
	double * getv_time();
	double * getv_amplitude();
	double getcharge();
	double gettime();
	
	//SETTERS
	void setcharge(double carica);
	void settime(double tempo);  
 
 private:

  double v_time[NSAMPLING];
  double v_amplitude[NSAMPLING];
  double charge;
  double time;
  
};

#endif
