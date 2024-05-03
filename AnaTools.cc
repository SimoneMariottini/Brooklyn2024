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
#include <iostream>
#include <cstdlib>
#include <string>
#include <string.h>
#include <TMinuit.h>
#include <math.h>
#include <TAxis.h>
#include <TString.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "AnaTools.h"

using namespace std;

//Constructor
AnaTools::AnaTools(TFile *f, Event *myEvent, double cf, double th){
  outfile_ = f;
  event_ = myEvent;
  cf_=cf;
  th_=th;

  nev_ = 0;

  return;
}

AnaTools::AnaTools(TFile *f){
  outfile_ = nullptr;
  event_ = nullptr;
  cf_=0;
  th_=0;
  nev_ = 0;

  TDirectory *dir;
  //aggiungere funzione bookable che controlla il tipo di constructor usato e non permette di fare bookings.
  
  //Check for charge histogram
  dir = f->GetDirectory("Hist_Charge");
  if(dir){
    dir->cd();

    for(int i = 0; i < NCHANNELS; i++){
     h_c_vector_[i] = (TH1D*)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));
    }

    h_c_tot_ = (TH1D*)gDirectory->Get("Hist_Total_Charge");

    bookings_[2] = 1;
  }

  //Check for ToF histogram
  dir = f->GetDirectory("Hist_ToF");
  if(dir){
    dir->cd();

    for(int i = 0; i < NCHANNELS; i++){
     h_TOF_cfm_[i] = (TH1D*)gDirectory->Get(Form("Hist_TOF_cfm_ch%i", i));
    }
    for(int i = 0; i < NCHANNELS; i++){
     h_TOF_ft_[i] = (TH1D*)gDirectory->Get(Form("Hist_TOF_ft_ch%i", i));
    }

    bookings_[3] = 1;
  }

  //Check for Persistence Maps
  dir = f->GetDirectory("Persistence_Map");
  if(dir){
    dir->cd();

    for(int i = 0; i < NCHANNELS; i++){
     persistence_vector_[i] = (TH2D*)gDirectory->Get(Form("Persistence_Map_Channel_%i", i));
    }

    bookings_[1] = 1;
  }

  gDirectory->cd("..");

  return;
}

//Destructor
AnaTools::~AnaTools(){
  return;
}


void AnaTools::BookWaveform(){
  bookings_[0] = 1;

  outfile_->cd("..");

  for (int i = 0; i < 20; i++){
    TString d_name = Form("Event_%d", i + 1);
    gDirectory->mkdir(d_name);
    gDirectory->cd(d_name);
    for(unsigned int k=0; k< NCHANNELS;k++){
      TString name = Form("Event_%d_Channel_%d", i + 1, k);
      TString title = Form("Event %d, Channel %d; time[ns]; Amplitude[V]", i + 1, k);
      h_wave_vector_[i][k] = new TH1D(name, title, 1024, 0, 1024*SAMPLINGPERIOD);
    }
    gDirectory->cd("..");	
  };

  gDirectory->cd("..");

  outfile_->cd();

  return;
}

void AnaTools::BookPersistence(){
  bookings_[1] = 1;

  outfile_->cd("..");

  gDirectory->mkdir("Persistence_Map");
  gDirectory->cd("Persistence_Map");

  for(unsigned int k=0; k<NCHANNELS;k++){
    TString name = Form("Persistence_Map_Channel_%d",k);
    TString title = Form("Persistence Map, channel %d; Time[ns]; Amplitude[V]", k);
    persistence_vector_[k] = new TH2D(name, title, 1024, 0, 1024*SAMPLINGPERIOD, 500, -0.3, 0.05);
  }
  gDirectory->cd("..");

  outfile_->cd();

  return;
}

void AnaTools::BookCharge(){
  bookings_[2] = 1;

  outfile_->cd("..");

  gDirectory->mkdir("Hist_Charge");
  gDirectory->cd("Hist_Charge");
  TString name = Form("Hist_Total_Charge");
  TString title = Form("Total Charge distribution; Charge[nC]; Counts[#]");
  h_c_tot_ = new TH1D(name, title, 500, -2., 2.);
  for(unsigned int k=0; k<NCHANNELS;k++){
    TString name = Form("Hist_Charge_Channel_%d",k);
    TString title = Form("Charge distribution, channel %d; Charge[nC]; Counts[#]", k);
    h_c_vector_[k] = new TH1D(name, title, 500, -0.01, 0.05);
  }
  gDirectory->cd("..");

  outfile_->cd();

  return;
}

void AnaTools::BookToF(){
  bookings_[3] = 1;

  outfile_->cd("..");

  gDirectory->mkdir("Hist_ToF");
  gDirectory->cd("Hist_ToF");
  for(int i=0;i<NCHANNELS;i++){
    //CONSTANT FRACTION HISTOGRAMS
    TString name = Form("Hist_TOF_cfm_ch%d",i);
    TString title = Form("Distribution of time of arrival using Constant Fraction Method, constant=%g, channel %d, not cut; Time of arrival[ns]; Counts [#]",cf_,i);
    h_TOF_cfm_[i] = new TH1D(name, title, 500, -10, 10);
    
    
    //FIXED THRESHOLD HISTOGRAMS
    name = Form("Hist_TOF_ft_ch%d",i);
    title = Form("Distribution of time of arrival using Fixed Threshold Method, threshold=%gV, channel %d, not cut; Time of arrival[ns]; Counts [#]",th_,i);
    h_TOF_ft_[i] = new TH1D(name, title, 500, -10, 10);
  }
  
  gDirectory->cd("..");

  outfile_->cd();

  return;
}

//Method for Data Analysis: Gets charges and light yields and puts them into histograms
void AnaTools::Process(int nevent){

  nev_= nevent;

  double tot_charge=0;
  double timeCFD[16];
  double timeFT[16];
  double charge[16];

  //CALCOLO CARICA E TEMPO DI ARRIVO
  for(unsigned int i =0; i < event_->getWaveforms().size(); i++){
    charge[i] = ComputeCharge(event_->getWaveforms().at(i));
    timeCFD[i] = ComputeTimeCFD(event_->getWaveforms().at(i),cf_);
    timeFT[i] = ComputeTimeFT(event_->getWaveforms().at(i),th_);
    tot_charge += charge[i];
  }

  event_->settot_charge(tot_charge);


  
  outfile_->cd();

  if (bookings_[2] == 1){ //Fills charge histograms
    h_c_tot_->Fill(tot_charge);
    for(unsigned int i =0; i < event_->getWaveforms().size(); i++){
      h_c_vector_[i]->Fill(charge[i]);
    }
  }


  if (bookings_[2] == 1){ //Fills ToF histograms
    double timeCFD_ref = 0.5*(timeCFD[0]+timeCFD[15]);
    double timeFT_ref = 0.5*(timeFT[0]+timeFT[15]);
    for(unsigned int i =0; i < event_->getWaveforms().size(); i++){
      double tofCFD = timeCFD[i]-timeCFD_ref;
      double tofFT = timeFT[i]-timeFT_ref;
      h_TOF_cfm_[i]->Fill(tofCFD);
      h_TOF_ft_[i]->Fill(tofFT);
    }
  }
  
  if (bookings_[0] == 1){ //Fills Waveform graphs
    for(unsigned int i =0; i < event_->getWaveforms().size(); i++){
      if(nev_>1 && nev_<=20){
        for(int isa=0;isa<1024;isa++){
	        h_wave_vector_[nev_-1][i]->SetBinContent(isa+1, ( event_->getWaveforms().at(i)->getv_amplitude())[isa]);
        }
      }
    }
  }

  if (bookings_[1] == 1){
    for(unsigned int i =0; i < event_->getWaveforms().size(); i++){
      for(int isa=0;isa<1024;isa++){
        persistence_vector_[i]->Fill((event_->getWaveforms().at(i)->getv_time())[isa], (event_->getWaveforms().at(i)->getv_amplitude())[isa]);
      }
    }
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



TFitResultPtr* AnaTools::FitCharge(){

  if (bookings_[2] == 1){
    static TFitResultPtr fitResults[NCHANNELS];
    int fitError = 0;
    TCanvas* c1 = new TCanvas();

    for (int i = 0; i < NCHANNELS; i++){
      TF1* gaus = new TF1("pedgaus", "[0]*TMath::Gaus(x,[1],[2])", -2, 2);
      TF1* langaus = new TF1("langaus", lanGausFun, -2, 2, 4);
      TF1* langausgaus = new TF1("langausgaus", lanGausPlusGausFun, -2, 2, 7);
      cout << "\r" <<"Fitting channel " << i << "." << flush;
      int nbins = h_c_vector_[i]->GetNbinsX();
      int midbin, infbin;
      double y1, y2, xmin, xmid;
      double par[4] = {0., 0., 0. ,0};
      double par2[7] = {0., 0., 0., 0., 0., 0., 0.};
      int j = 1;

      //Looks for the bin that isn't almost empty
      while(j <= nbins){
        if(h_c_vector_[i]->GetBinContent(j) > 10){
          xmin = h_c_vector_[i]->GetBinCenter(j - 5);
          infbin = j-5;
          break;
        }
        j++;
      }
      
      //Looks for the first maximum bin

      while(j <= nbins){
        y1 = h_c_vector_[i]->GetBinContent(j);
        y2 = h_c_vector_[i]->GetBinContent(j + 1);
        if(y2 < y1){
          break;
        }
        j++;
      }

      //Looks for the local minimum bin after the first maximum 

      while(j <= nbins){
        y1 = h_c_vector_[i]->GetBinContent(j);
        y2 = h_c_vector_[i]->GetBinContent(j + 1);
        if(y2 > y1){
          break;
        }
        j++;
      }
      //save local minimum bin number and position
      midbin = j; 
      xmid = h_c_vector_[i]->GetBinCenter(midbin);
      h_c_vector_[i]->GetXaxis()->SetRange(1, midbin);
      
      gaus->SetParameters(h_c_vector_[i]->Integral(), h_c_vector_[i]->GetMean(), h_c_vector_[i]->GetStdDev());

      gaus->SetParNames("Constant", "Mean", "Sigma");

      h_c_vector_[i]->Fit("pedgaus", "Q", "same", xmin, xmid); //Fitting the pedestal gaussian

      //Saving the pedestal gaussian fit parameters
      par2[4] = gaus->GetParameter("Constant");
      par2[5] = gaus->GetParameter("Mean");
      par2[6] = gaus->GetParameter("Sigma");

      h_c_vector_[i]->GetXaxis()->SetRange(0, 0); //reset range

      //Temporarely subtract the pedestal gaussian to better fit the langaus

      vector<int>* reallocBin = new vector<int>;
      vector<double>* reallocBinContent = new vector<double>;


      for(int k = infbin; k < midbin + 20; k++){
        double adjValue = h_c_vector_[i]->GetBinContent(k) - gaus->Eval(h_c_vector_[i]->GetBinCenter(k));
        reallocBin->push_back(k);
        reallocBinContent->push_back(h_c_vector_[i]->GetBinContent(k));
        if(adjValue >= 0){
          h_c_vector_[i]->SetBinContent(k, adjValue);
        }
        else{
          h_c_vector_[i]->SetBinContent(k, 0.);
        }
      }

      //Evaluating the best starting parameters for the langaus

      par[2] = h_c_vector_[i]->Integral(); //Normalization factor

      j = h_c_vector_[i]->GetMaximumBin();
      par[1] = h_c_vector_[i]->GetBinCenter(j); //Most Probable value

      par[0] = h_c_vector_[i]->GetStdDev()*0.1; //Width

      par[3] = 0.01; //Gaussian resolution

      //Fit of langaus function
      langaus->SetNpx(1000);
      langaus->SetParameters(par);
      langaus->SetParLimits(3, 0.001, 10);
      langaus->SetParNames("Width", "MP", "Area", "GSigma");

      h_c_vector_[i]->Fit(langaus, "Q", "same", -0.01, 0.05);

      //Saving the langaus fit parameters
      par2[0] = langaus->GetParameter("Width");
      par2[1] = langaus->GetParameter("MP");
      par2[2] = langaus->GetParameter("Area");
      par2[3] = langaus->GetParameter("GSigma");

      //Returning the histogram to normal
      for(int k = 0; k < reallocBin->size(); k++){
        h_c_vector_[i]->SetBinContent((*reallocBin)[k], (*reallocBinContent)[k]);
      }
      delete reallocBin;
      delete reallocBinContent;
      
      //Fit of total distribution using the previous fit parameters (pedestal gaus + langaus)
      langausgaus->SetNpx(1000);
      langausgaus->SetParameters(par2);
      langausgaus->SetParLimits(3, 0.001, 10);
      langausgaus->SetParNames("Width", "MP", "Area", "G Sigma", "Ped Constant", "Ped Mean", "Ped Sigma");

      fitResults[i] = h_c_vector_[i]->Fit(langausgaus, "QS", "same", -0.01, 0.1);
      
      if((int)fitResults[i] != 0){
        fitError++;
        cout << "Fit of channel " << i << " was unsuccessful! Exit status " << (int)fitResults[i] << endl;
      }

      delete gaus;
      delete langaus;
      delete langausgaus;
    };
    delete c1;
    cout << "\r" << NCHANNELS - fitError << "/" << NCHANNELS << " fits were successful" << endl;
    return fitResults; //returns the number of unsuccessfull fits
  }
  else{
    cout << "Can't fit charge histograms." << endl;
    return nullptr; 
  }
  
}

double AnaTools::lanGausFun(double *x, double *par) { //copied from laungaus.C example in Root documentation
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double relsigma = 0.000000001;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      //Relative value of Gaussian sigma
      //relsigma = par[3]*0.001; //questo funziona ma praticamente la gaussiana è da trascurare
      if(x[0] >= 0){
        relsigma = par[3]*TMath::Sqrt(TMath::Abs(x[0]));
       
 
      // Range of convolution integral
      xlow = x[0] - sc * relsigma;
      xupp = x[0] + sc * relsigma;
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,relsigma);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,relsigma);
      }
      return (par[2] * step * sum * invsq2pi / relsigma);
      }

      else{
        return (par[2]*TMath::Landau(x[0],mpc,par[0])/par[0]);
      }
 
      
}

double AnaTools::lanGausPlusGausFun(double *x, double *par){
  double lanpar[4] = {par[0], par[1], par[2], par[3]};
  return lanGausFun(x, lanpar) + par[4]*TMath::Gaus(x[0],par[5],par[6]);
}

double* AnaTools::EvaluateEfficiency(){
  int goodFits = 0;
  static double result[NCHANNELS*2] = {0};
  TFitResultPtr* fitResults = FitCharge();

  for(int i = 0; i < NCHANNELS; i++){
    if((int)fitResults[i] == 0) goodFits++;
  }
  
  if(goodFits != NCHANNELS){
    cout << "Not all fits were successful! Can't evaluate channel efficiency." << endl;
    return nullptr;
  }

  for(int i = 0; i < NCHANNELS; i++){
    double pedArea, langausArea, sigmaPed, sigmaLangaus;
    double* param;
    //double* covMatrix;
    //retrieve the the different fit functions
    TF1* langausgaus = h_c_vector_[i]->GetFunction("langausgaus");
    TF1* pedgaus = new TF1("pedgaus", "[0]*TMath::Gaus(x,[1],[2])", -2, 2);
    TF1* langaus = new TF1("langaus", lanGausFun, -2, 2, 4);

    //Get fit parameters
    
    param = langausgaus->GetParameters();
    //covMatrix = fitResults[i]->GetCovarianceMatrix().GetMatrixArray();

    TMatrixDSym covMatrix = fitResults[i]->GetCovarianceMatrix();

    //Set the final fit parameters
    langaus->SetParameters(param);

    pedgaus->SetParameters(param + 4);

    //Calculate the area under the two functions
    pedArea = pedgaus->Integral(-0.01, 0.05);
    sigmaPed = pedgaus->IntegralError(-0.01, 0.05, param + 4, covMatrix.GetSub(4, 6, 4, 6).GetMatrixArray());

    langausArea = langaus->Integral(-0.01, 0.05);
    sigmaLangaus = langaus->IntegralError(-0.01, 0.05, param, covMatrix.GetSub(0, 3, 0, 3).GetMatrixArray());

    result[i] = langausArea/(langausArea + pedArea);
    result[i + NCHANNELS] = 1/TMath::Power((langausArea + pedArea), 2)*TMath::Sqrt(TMath::Power(pedArea, 2)*TMath::Power(sigmaLangaus, 2) + TMath::Power(langausArea, 2)*TMath::Power(sigmaPed, 2));

    delete pedgaus;
    delete langaus;
  }

  return result;
}


   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
   	
 
