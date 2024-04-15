#ifndef EVENT_H
#define EVENT_H

#include <stdint.h>
#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include "Waveform.h"

using namespace std;

class Event{
  
 public:
  
  Event();
  virtual ~Event();

  void ReadEvent(string inname,unsigned long int* pos);
  void Clear();

 //GETTERS
  vector<Waveform*> getWaveforms(); 
  double GetTot_charge();
  double gettime();
  
  //SETTERS
  void settot_charge(double caricaTotale);
  void settime(double tempo);
  
 private:

  vector<Waveform *>myWaveforms;
  double tot_charge;
  double time;
  
};

#endif
