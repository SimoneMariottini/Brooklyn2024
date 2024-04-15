#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "Waveform.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <TMinuit.h>
#include <math.h>


using namespace std;

//COSTRUTTORE E DISTRUTTORE
Waveform::Waveform(double *vt, double *va){
  charge =-1000;
  time = -1000;
  //set time and amplitudes array with the values read from file
  memcpy(v_time, vt, sizeof(v_time));   //memcpy mette vt e va nella locazione di memoria di v_time e v_amplitude 
  memcpy(v_amplitude, va, sizeof(v_amplitude));
}
Waveform::~Waveform(){};

//GETTERS
double* Waveform::getv_time(){
	return v_time;
}
double* Waveform::getv_amplitude(){
	return v_amplitude;
}

double Waveform::getcharge(){
	return charge;
}
double Waveform::gettime(){
	return time;
}
//SETTERS
void Waveform::setcharge(double carica){
	charge=carica;	
}
void Waveform::settime(double tempo){
	time=tempo;
}
