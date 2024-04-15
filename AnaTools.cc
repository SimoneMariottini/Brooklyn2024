#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "AnaTools.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <string.h>
#include <TMinuit.h>
#include <math.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TColor.h>
#include <TCanvas.h>

using namespace std;

//Constructor
AnaTools::AnaTools(TFile *f, Event *myEvent, double cf_, double th_){ 

  outfile = f;
  event = myEvent;
  cf=cf_;
  th=th_;

  nev = 0;
  
}




//Destructor
AnaTools::~AnaTools(){
  return;

}





void AnaTools::BookingHistograms(){

  outfile->cd("..");

  //PERSISTENCE MAP
  gDirectory->mkdir("Persistence_Map");
  gDirectory->cd("Persistence_Map");

  for(unsigned int k=0; k<NCHANNELS;k++){
    TString name = Form("Persistence_Map_Channel_%d",k);
    TString title = Form("Persistence Map, channel %d; Time[]; Amplitude[]", k);
    persistence_vector[k] = new TH2D(name, title, 1024, 0, 1024, 500, -0.3, 0.05);
  }
  gDirectory->cd("..");

  //CHARGE HISTOGRAMS
  gDirectory->mkdir("Hist_Charge");
  gDirectory->cd("Hist_Charge");
  TString name = Form("Hist_total_charge");
  TString title = Form("Total Charge distribution; Charge[nC]; Counts[#]");
  hctot = new TH1D(name, title, 500, -2., 2.);
  for(unsigned int k=1; k<=NCHANNELS;k++){
    TString name = Form("Hist_Channel_%d",k-1);
    TString title = Form("Charge distribution, channel %d; Charge[nC]; Counts[#]", k-1);
    hc_vector[k-1] = new TH1D(name, title, 500, -0.1, 0.1);
  }
  gDirectory->cd("..");

    
   
  //TOF HISTOGRAMS
  gDirectory->mkdir("Hist_ToF");
  gDirectory->cd("Hist_ToF");
  for(int i=1;i<=NCHANNELS;i++){
    //CONSTANT FRACTION HISTOGRAMS
    
    TString name = Form("Hist_TOF_cfm_ch%d",i);
    TString title = Form("Distribution of time of arrival using Constant Fraction Method, constant=%g, channel %d, not cut; Time of arrival[ns]; Counts [#]",cf,i);
    hTOF_cfm[i-1] = new TH1D(name, title, 500, -10, 10);
    
    
    //FIXED THRESHOLD HISTOGRAMS
    name = Form("Hist_TOF_ft_ch%d",i);
    title = Form("Distribution of time of arrival using Fixed Threshold Method, threshold=%gV, channel %d, not cut; Time of arrival[ns]; Counts [#]",th,i);
    hTOF_ft[i-1] = new TH1D(name, title, 500, -10, 10);
  }
  
  gDirectory->cd("..");
 
  
  for(int i=1; i<=20; i++){       
    //WAVEFORMS HISTOGRAMS                
    char dir[10]="Evento";
    char *newEvent=strcat(dir, to_string(i).c_str());        
    gDirectory->mkdir(&newEvent[0]);
    gDirectory->cd(&newEvent[0]);
    for(unsigned int k=1; k<=NCHANNELS;k++){
      TString name = Form("Event_%d_Channel_%d", i,k-1);
      TString title = Form("Event %d, Channel %d; time[ns]; Amplitude[V]", i,k-1);
      hist_vector[i-1][k-1] = new TH1D(name, title, 1024, 0, 1024*SAMPLINGPERIOD);
    }
    gDirectory->cd("..");	
  }

  
 
  gDirectory->cd("..");
	
	
}

//Method for Data Analysis: Gets charges and light yields and puts them into histograms
void AnaTools::Process(int nevent){

  nev= nevent;



  double tot_charge=0;
  double timeCFD[16];
  double timeFT[16];
  double charge[16];

  //CALCOLO CARICA E TEMPO DI ARRIVO
  for(unsigned int i =0; i < event->getWaveforms().size(); i++){
    charge[i] = ComputeCharge(event->getWaveforms().at(i));
    timeCFD[i] = ComputeTimeCFD(event->getWaveforms().at(i),cf);
    timeFT[i] = ComputeTimeFT(event->getWaveforms().at(i),th);
    tot_charge += charge[i];
  }

  event->settot_charge(tot_charge);


  
  outfile->cd();

  //CHARGE HISTOGRAMS FILLING
  hctot->Fill(tot_charge);
  for(unsigned int i =0; i < event->getWaveforms().size(); i++){
    hc_vector[i]->Fill(charge[i]);
  }



  //TOF HISTOGRAMS FILLING
  double timeCFD_ref = 0.5*(timeCFD[0]+timeCFD[15]);
  double timeFT_ref = 0.5*(timeFT[0]+timeFT[15]);
  for(unsigned int i =0; i < event->getWaveforms().size(); i++){
    double tofCFD = timeCFD[i]-timeCFD_ref;
    double tofFT = timeFT[i]-timeFT_ref;
    hTOF_cfm[i]->Fill(tofCFD);
    hTOF_ft[i]->Fill(tofFT);
  }

  
  //WAVEFORM HISTOGRAMS FILLING
  for(unsigned int i =0; i < event->getWaveforms().size(); i++){
    if(nev>=1 && nev<=20){
      for(int isa=0;isa<1024;isa++){
	hist_vector[nev-1][i]->SetBinContent(isa+1, ( event->getWaveforms().at(i)->getv_amplitude())[isa]);
      }
    }
  }

  //PERSISTENCE HISTOGRAM FILLING
  //if(nev>=1 && nev<=100){
  for(unsigned int i =0; i < event->getWaveforms().size(); i++){
    for(int isa=0;isa<1024;isa++){
      persistence_vector[i]->Fill((double)isa, (event->getWaveforms().at(i)->getv_amplitude())[isa]);
    }
  //}
  }

  return;
}






double AnaTools::ComputeCharge(Waveform *w){

  double *tmp_amp = w->getv_amplitude();
  double *tmp_time = w->getv_time();

  double charge = 0.;
  double prod;
  double chargetot=0;

  double deltat=0.; //trapezi
  double b1=0, b2=0;
  double a1, a2;
  double m,q,t1,t2,thalf, dq;

  for (int j = 0; j<1023; j++){
    a1 = tmp_amp[j];
    a2 = tmp_amp[j+1];
    t1 = tmp_time[j];
    t2 = tmp_time[j+1];
    prod=a1*a2;
    if(prod<0){
      m = (a2-a1)/(t2-t1);
      q = a2-m*t2;
      thalf = -q/m;
      dq=0.5*a1*(thalf-t1)+0.5*a2*(t2-thalf);
    }else{
      dq = 0.5*(a1+a2)*(t2-t1);
    }

    
    charge+=dq;
  }

  //convert to nC (divede 50ohm)
  charge*=-1;
  charge/=50;
  

  
  return charge;
	
}




double AnaTools::ComputeTimeCFD(Waveform *w, double fraction){

  double *tmp_amp = w->getv_amplitude();
  double *tmp_time = w->getv_time();
  
  //evaluate the baseline (media dei primi 50 campionamenti)
  double baseline=0;
  for(int i=0;i<20;i++)baseline+=tmp_amp[i]/20;


  //evaluate the minimum amplitude (attenzione e' negativa)
  double minamp=0;
  int i_ampmin=0;
  for(int i=0;i<1024;i++){
    if(tmp_amp[i]<minamp){
      minamp=tmp_amp[i];
      i_ampmin=i;
    }
  }
  //prendo il modulo
  double maxamp = fabs(minamp);

  
  //evaluate the absolute threshold
  double AbsoluteThreshold=-fraction*maxamp+baseline;

  Int_t i_thr = i_ampmin;

  double t_arr=-1000;
  bool foundthreshold = false;
  while(!foundthreshold && i_thr<1023 && i_thr>1){
    double a1 = tmp_amp[i_thr];
    double a2 = tmp_amp[i_thr+1];
    double t1 = tmp_time[i_thr];
    double t2 = tmp_time[i_thr+1];
    double m = (a2-a1)/(t2-t1);
    double q = a1 - m*t1;
    t_arr = (AbsoluteThreshold-q)/m;
    if(AbsoluteThreshold>a2 && AbsoluteThreshold<=a1)foundthreshold =true;
    i_thr--;
  }


  
  return t_arr;
  
	
}




double AnaTools::ComputeTimeFT(Waveform *w, double threshold){

   double *tmp_amp = w->getv_amplitude();
  double *tmp_time = w->getv_time();
  
  //evaluate the baseline (media dei primi 50 campionamenti)
  double baseline=0;
  for(int i=0;i<20;i++)baseline+=tmp_amp[i]/20;


  //evaluate the minimum amplitude (attenzione e' negativa)
  double minamp=0;
  int i_ampmin=0;
  for(int i=0;i<1024;i++){
    if(tmp_amp[i]<minamp){
      minamp=tmp_amp[i];
      i_ampmin=i;
    }
  }
  //prendo il modulo
  double maxamp = fabs(minamp);

  
  //evaluate the absolute threshold
  double AbsoluteThreshold=-threshold+baseline;

  Int_t i_thr = i_ampmin;

  double t_arr=-1000;
  bool foundthreshold = false;
  while(!foundthreshold && i_thr<1023 && i_thr>1){
    double a1 = tmp_amp[i_thr];
    double a2 = tmp_amp[i_thr+1];
    double t1 = tmp_time[i_thr];
    double t2 = tmp_time[i_thr+1];
    double m = (a2-a1)/(t2-t1);
    double q = a1 - m*t1;
    t_arr = (AbsoluteThreshold-q)/m;
    if(AbsoluteThreshold>a2 && AbsoluteThreshold<=a1)foundthreshold =true;
    i_thr--;
  }
  
  return t_arr;
	

	
}












void AnaTools::LoadPedestal(string inname){
  // ifstream infile; //File "inname.dat" stream
  // string fileline; //A line of the File
  // int ch=0;
  // infile.open(inname);// file containing numbers in 5 columns
  // getline(infile,fileline);
  // if(infile.fail()) // checks to see if file opended 
  //   { 
  //     cout << "error= input file of pedestal not found." << endl; 
  //   } 
  // while(!infile.eof()) // reads file to end of file
  //   { 
 						
  //     getline(infile,fileline);
  //     istringstream ss(fileline);
  //     ss >> ped_mean[ch];
  //     ss >> ped_mean[ch];
  //     ss >> ped_sigma[ch];
  //     ss >> mean2[ch];
  //     ss >> sigma2[ch];
  //     ch++;
  //   }
  // //EVALUATING GAIN    
  // for(int ch=0; ch<NCHANNELS; ch++){
  //   gain_[ch]=mean2[ch]/ECHARGE;
  // }
  // infile.close();

  /*for(int i=0;i<16;i++)
    { 
    cout << "\tped mean= "<< ped_mean[i] << "\tped sigma= " << ped_sigma[i] << endl;
   	
   	 	
    }*/
  
}   	

void AnaTools::Draw(){
  gStyle->SetPalette(kRainBow);
  for(unsigned int i =0; i < event->getWaveforms().size(); i++){
    auto c1 = new TCanvas("c1","c1",600,400);
    //persistence_vector[i]->SetContour(255);
    persistence_vector[i]->Draw("CONT2"); 
    delete c1;
  }
}


   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
 
