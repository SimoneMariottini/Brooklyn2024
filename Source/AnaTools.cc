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

#include "../Header/AnaTools.h"

using namespace std;

// Constructor
AnaTools::AnaTools(TFile *f, Event *myEvent, double cf, double th)
{
  outfile_ = f;
  event_ = myEvent;
  cf_ = cf;
  th_ = th;

  nev_ = 0;

  return;
}

AnaTools::AnaTools(TFile *f)
{
  outfile_ = nullptr;
  event_ = nullptr;
  cf_ = 0;
  th_ = 0;
  nev_ = 0;

  TDirectory *dir;
  // aggiungere funzione bookable che controlla il tipo di constructor usato e non permette di fare bookings.

  // Check for charge histogram
  dir = f->GetDirectory("Hist_Charge");
  if (dir)
  {
    dir->cd();

    for (int i = 0; i < NCHANNELS; i++)
    {
      h_c_vector_[i] = (TH1D *)gDirectory->Get(Form("Hist_Charge_Channel_%i", i));
    }

    h_c_tot_ = (TH1D *)gDirectory->Get("Hist_Total_Charge");

    bookings_[2] = 1;
  }

  // Check for ToF histogram
  dir = f->GetDirectory("Hist_ToF");
  if (dir)
  {
    dir->cd();

    for (int i = 0; i < NCHANNELS; i++)
    {
      h_TOF_cfm_[i] = (TH1D *)gDirectory->Get(Form("Hist_TOF_cfm_ch%i", i));
    }
    for (int i = 0; i < NCHANNELS; i++)
    {
      h_TOF_ft_[i] = (TH1D *)gDirectory->Get(Form("Hist_TOF_ft_ch%i", i));
    }

    bookings_[3] = 1;
  }

  // Check for Persistence Maps
  dir = f->GetDirectory("Persistence_Map");
  if (dir)
  {
    dir->cd();

    for (int i = 0; i < NCHANNELS; i++)
    {
      persistence_vector_[i] = (TH2D *)gDirectory->Get(Form("Persistence_Map_Channel_%i", i));
    }

    bookings_[1] = 1;
  }

  gDirectory->cd("..");

  return;
}

AnaTools::AnaTools(TFile *f, Event *myEvent, double cf, double th, TString infoFile) : AnaTools::AnaTools(f, myEvent, cf, th)
{
  infoFile_ = infoFile;
  LoadInfo(infoFile_);
  return;
}

AnaTools::AnaTools(TFile *f, TString infoFile) : AnaTools::AnaTools(f)
{
  infoFile_ = infoFile;
  LoadInfo(infoFile_);
  return;
}

// Destructor
AnaTools::~AnaTools()
{
  return;
}

void AnaTools::BookWaveform()
{
  bookings_[0] = 1;

  outfile_->cd("..");

  for (int i = 0; i < 20; i++)
  {
    TString d_name = Form("Event_%d", i + 1);
    gDirectory->mkdir(d_name);
    gDirectory->cd(d_name);
    for (unsigned int k = 0; k < NCHANNELS; k++)
    {
      TString name = Form("Event_%d_Channel_%d", i + 1, k);
      TString title = Form("Event %d, Channel %d; time[ns]; Amplitude[V]", i + 1, k);
      h_wave_vector_[i][k] = new TH1D(name, title, 1024, 0, 1024 * SAMPLINGPERIOD);
    }
    gDirectory->cd("..");
  };

  gDirectory->cd("..");

  outfile_->cd();

  return;
}

void AnaTools::BookPersistence()
{
  bookings_[1] = 1;

  outfile_->cd("..");

  gDirectory->mkdir("Persistence_Map");
  gDirectory->cd("Persistence_Map");

  for (unsigned int k = 0; k < NCHANNELS; k++)
  {
    TString name = Form("Persistence_Map_Channel_%d", k);
    TString title = Form("Persistence Map, channel %d; Time[ns]; Amplitude[V]", k);
    persistence_vector_[k] = new TH2D(name, title, 1024, 0, 1024 * SAMPLINGPERIOD, 500, -0.3, 0.05);
  }
  gDirectory->cd("..");

  outfile_->cd();

  return;
}

void AnaTools::BookCharge()
{
  bookings_[2] = 1;

  outfile_->cd("..");

  gDirectory->mkdir("Hist_Charge");
  gDirectory->cd("Hist_Charge");
  TString name = Form("Hist_Total_Charge");
  TString title = Form("Total Charge distribution; Charge[nC]; Counts[#]");
  h_c_tot_ = new TH1D(name, title, 500, -2., 2.);
  for (unsigned int k = 0; k < NCHANNELS; k++)
  {
    TString name = Form("Hist_Charge_Channel_%d", k);
    TString title = Form("Charge distribution, channel %d; Charge[nC]; Counts[#]", k);
    h_c_vector_[k] = new TH1D(name, title, 500, h_charge_range_[0], h_charge_range_[1]);
  }
  gDirectory->cd("..");

  outfile_->cd();

  return;
}

void AnaTools::BookCharge(const double& xmin, const double& xmax){ //Book charge setting custom histogram range
  this->SetChargeRange(xmin, xmax);
  this->BookCharge();
  return;
}

void AnaTools::BookToF(){
  bookings_[3] = 1;

  outfile_->cd("..");

  gDirectory->mkdir("Hist_ToF");
  gDirectory->cd("Hist_ToF");
  for (int i = 0; i < NCHANNELS; i++)
  {
    // CONSTANT FRACTION HISTOGRAMS
    TString name = Form("Hist_TOF_cfm_ch%d", i);
    TString title = Form("Distribution of time of arrival using Constant Fraction Method, constant=%g, channel %d, not cut; Time of arrival[ns]; Counts [#]", cf_, i);
    h_TOF_cfm_[i] = new TH1D(name, title, 500, -10, 10);

    // FIXED THRESHOLD HISTOGRAMS
    name = Form("Hist_TOF_ft_ch%d", i);
    title = Form("Distribution of time of arrival using Fixed Threshold Method, threshold=%gV, channel %d, not cut; Time of arrival[ns]; Counts [#]", th_, i);
    h_TOF_ft_[i] = new TH1D(name, title, 500, -10, 10);
  }

  gDirectory->cd("..");

  outfile_->cd();

  return;
}

void AnaTools::BookTime()
{

  bookings_[4] = 1;

  outfile_->cd("..");

  gDirectory->mkdir("Hist_Time");
  gDirectory->cd("Hist_Time");

  TString name = "Hist_Time";
  TString title = "Distribution of time of arrival";

  for(int j = 1; j <= 2; j++){
    gDirectory->mkdir(Form("Hist_Time_Method_%i", j));
    gDirectory->cd(Form("Hist_Time_Method_%i", j));

    for(int i = 0; i < NCHANNELS; i++){

      TString name = Form("Hist_time_Method_%i_Channel_%i", j, i);
      TString title = Form("Distribution of time of arrival, method %i, channel %i", j, i);
      
      h_time_vector_[j - 1][i] = new TH1D(name, title, 500, -100, 100);
    }
    gDirectory->cd("..");
  }

  h_time_ = new TH1D(name, title, 500, -100, 100);

  gDirectory->cd("..");

  outfile_->cd();

  return;
}

void AnaTools::Process(int nevent)
{

  nev_ = nevent;

  double tot_charge = 0;
  double timeCFD[16];
  double timeFT[16];
  double charge[16];

  // CALCOLO CARICA E TEMPO DI ARRIVO
  for (unsigned int i = 0; i < event_->getWaveforms().size(); i++)
  {
    charge[i] = ComputeCharge(event_->getWaveforms().at(i));
    timeCFD[i] = ComputeTimeCFD(event_->getWaveforms().at(i), cf_);
    timeFT[i] = ComputeTimeFT(event_->getWaveforms().at(i), th_);
    tot_charge += charge[i];
  }

  event_->settot_charge(tot_charge);

  outfile_->cd();

  if (bookings_[2] == 1)
  { // Fills charge histograms
    h_c_tot_->Fill(tot_charge);
    for (unsigned int i = 0; i < event_->getWaveforms().size(); i++)
    {
      h_c_vector_[i]->Fill(charge[i]);
    }
  }

  if (bookings_[3] == 1)
  { // Fills ToF histograms
    double timeCFD_ref = 0.5 * (timeCFD[0] + timeCFD[15]);
    double timeFT_ref = 0.5 * (timeFT[0] + timeFT[15]);
    for (unsigned int i = 0; i < event_->getWaveforms().size(); i++)
    {
      double tofCFD = timeCFD[i] - timeCFD_ref;
      double tofFT = timeFT[i] - timeFT_ref;
      h_TOF_cfm_[i]->Fill(tofCFD);
      h_TOF_ft_[i]->Fill(tofFT);
    }
  }

  if (bookings_[4] == 1){

    int novercutoff = 0;
    double TimeInternal = 0;
    double TimeExternal = 0.5 * (timeCFD[0] + timeCFD[NCHANNELS - 1]);
    double PartialInternalTime, totalTime;

    for (int i = 1; i < NCHANNELS - 1; i++)
    {
      if (charge[i] > cutoff_[i])
      {
        TimeInternal += timeCFD[i];
        novercutoff++;

        h_time_vector_[0][i]->Fill(timeCFD[i] - TimeExternal);
      }
    }
    totalTime = TimeInternal + TimeExternal * 2;
    for (int i = 1; i < NCHANNELS - 1; i++)
    {
      if (charge[i] > cutoff_[i])
      {
        PartialInternalTime = totalTime - timeCFD[i];
        h_time_vector_[1][i]->Fill(PartialInternalTime / (novercutoff + 1) - timeCFD[i]);
      }
    }

    TimeInternal /= novercutoff;
    h_time_->Fill(TimeExternal - TimeInternal);
  }

  if (bookings_[0] == 1){ // Fills Waveform graphs
    for (unsigned int i = 0; i < event_->getWaveforms().size(); i++)
    {
      if (nev_ > 1 && nev_ <= 20)
      {
        for (int isa = 0; isa < 1024; isa++)
        {
          h_wave_vector_[nev_ - 1][i]->SetBinContent(isa + 1, (event_->getWaveforms().at(i)->getv_amplitude())[isa]);
        }
      }
    }
  }

  if (bookings_[1] == 1)
  {
    for (unsigned int i = 0; i < event_->getWaveforms().size(); i++)
    {
      for (int isa = 0; isa < 1024; isa++)
      {
        persistence_vector_[i]->Fill((event_->getWaveforms().at(i)->getv_time())[isa], (event_->getWaveforms().at(i)->getv_amplitude())[isa]);
      }
    }
  }

  return;
}

double AnaTools::ComputeCharge(Waveform *w)
{

  double *tmp_amp = w->getv_amplitude();
  double *tmp_time = w->getv_time();

  double charge = 0.;
  double prod;
  double chargetot = 0;

  double deltat = 0.; // trapezi
  double b1 = 0, b2 = 0;
  double a1, a2;
  double m, q, t1, t2, thalf, dq;

  for (int j = 0; j < 1023; j++)
  {
    a1 = tmp_amp[j];
    a2 = tmp_amp[j + 1];
    t1 = tmp_time[j];
    t2 = tmp_time[j + 1];
    prod = a1 * a2;
    if (prod < 0)
    {
      m = (a2 - a1) / (t2 - t1);
      q = a2 - m * t2;
      thalf = -q / m;
      dq = 0.5 * a1 * (thalf - t1) + 0.5 * a2 * (t2 - thalf);
    }
    else
    {
      dq = 0.5 * (a1 + a2) * (t2 - t1);
    }

    charge += dq;
  }

  // convert to nC (divede 50ohm)
  charge *= -1;
  charge /= 50;

  return charge;
}

double AnaTools::ComputeTimeCFD(Waveform *w, double fraction)
{

  double *tmp_amp = w->getv_amplitude();
  double *tmp_time = w->getv_time();

  // evaluate the baseline (media dei primi 50 campionamenti)
  double baseline = 0;
  for (int i = 0; i < 20; i++)
    baseline += tmp_amp[i] / 20;

  // evaluate the minimum amplitude (attenzione e' negativa)
  double minamp = 0;
  int i_ampmin = 0;
  for (int i = 0; i < 1024; i++)
  {
    if (tmp_amp[i] < minamp)
    {
      minamp = tmp_amp[i];
      i_ampmin = i;
    }
  }
  // prendo il modulo
  double maxamp = fabs(minamp);

  // evaluate the absolute threshold
  double AbsoluteThreshold = -fraction * maxamp + baseline;

  Int_t i_thr = i_ampmin;

  double t_arr = -1000;
  bool foundthreshold = false;
  while (!foundthreshold && i_thr < 1023 && i_thr > 1)
  {
    double a1 = tmp_amp[i_thr];
    double a2 = tmp_amp[i_thr + 1];
    double t1 = tmp_time[i_thr];
    double t2 = tmp_time[i_thr + 1];
    double m = (a2 - a1) / (t2 - t1);
    double q = a1 - m * t1;
    t_arr = (AbsoluteThreshold - q) / m;
    if (AbsoluteThreshold > a2 && AbsoluteThreshold <= a1)
      foundthreshold = true;
    i_thr--;
  }

  return t_arr;
}

double AnaTools::ComputeTimeFT(Waveform *w, double threshold)
{

  double *tmp_amp = w->getv_amplitude();
  double *tmp_time = w->getv_time();

  // evaluate the baseline (media dei primi 50 campionamenti)
  double baseline = 0;
  for (int i = 0; i < 20; i++)
    baseline += tmp_amp[i] / 20;

  // evaluate the minimum amplitude (attenzione e' negativa)
  double minamp = 0;
  int i_ampmin = 0;
  for (int i = 0; i < 1024; i++)
  {
    if (tmp_amp[i] < minamp)
    {
      minamp = tmp_amp[i];
      i_ampmin = i;
    }
  }
  // prendo il modulo
  double maxamp = fabs(minamp);

  // evaluate the absolute threshold
  double AbsoluteThreshold = -threshold + baseline;

  Int_t i_thr = i_ampmin;

  double t_arr = -1000;
  bool foundthreshold = false;
  while (!foundthreshold && i_thr < 1023 && i_thr > 1)
  {
    double a1 = tmp_amp[i_thr];
    double a2 = tmp_amp[i_thr + 1];
    double t1 = tmp_time[i_thr];
    double t2 = tmp_time[i_thr + 1];
    double m = (a2 - a1) / (t2 - t1);
    double q = a1 - m * t1;
    t_arr = (AbsoluteThreshold - q) / m;
    if (AbsoluteThreshold > a2 && AbsoluteThreshold <= a1)
      foundthreshold = true;
    i_thr--;
  }

  return t_arr;
}

TFitResultPtr *AnaTools::FitCharge()
{

  if (bookings_[2] == 1)
  {
    static TFitResultPtr fitResults[NCHANNELS];
    int fitError = 0;
    TCanvas *c1 = new TCanvas();

    for (int i = 0; i < NCHANNELS; i++)
    {
      TF1 *gaus = new TF1("pedgaus", "[0]*TMath::Gaus(x,[1],[2])", -2, 2);
      TF1 *langaus = new TF1("langaus", lanGausFun, -2, 2, 4);
      TF1 *langausgaus = new TF1("langausgaus", lanGausPlusGausFun, -2, 2, 7);
      cout << "\r" << "Fitting channel " << i << "." << flush;
      int nbins = h_c_vector_[i]->GetNbinsX();
      int midbin, infbin;
      double y1, y2, xmin, xmid;
      double par[4] = {0., 0., 0., 0};
      double par2[7] = {0., 0., 0., 0., 0., 0., 0.};
      int j = 1;

      // Looks for the bin that isn't almost empty
      while (j <= nbins)
      {
        if (h_c_vector_[i]->GetBinContent(j) > 10)
        {
          xmin = h_c_vector_[i]->GetBinCenter(j - 5);
          infbin = j - 5;
          break;
        }
        j++;
      }

      // Looks for the first maximum bin

      while (j <= nbins)
      {
        y1 = h_c_vector_[i]->GetBinContent(j);
        y2 = h_c_vector_[i]->GetBinContent(j + 1);
        if (y2 < y1)
        {
          break;
        }
        j++;
      }

      // Looks for the local minimum bin after the first maximum

      while (j <= nbins)
      {
        y1 = h_c_vector_[i]->GetBinContent(j);
        y2 = h_c_vector_[i]->GetBinContent(j + 1);
        if (y2 > y1)
        {
          break;
        }
        j++;
      }
      // save local minimum bin number and position
      midbin = j;
      xmid = h_c_vector_[i]->GetBinCenter(midbin);
      h_c_vector_[i]->GetXaxis()->SetRange(1, midbin);

      gaus->SetParameters(h_c_vector_[i]->Integral(), h_c_vector_[i]->GetMean(), h_c_vector_[i]->GetStdDev());

      gaus->SetParNames("Constant", "Mean", "Sigma");

      h_c_vector_[i]->Fit("pedgaus", "Q", "same", xmin, xmid); // Fitting the pedestal gaussian

      // Saving the pedestal gaussian fit parameters
      par2[4] = gaus->GetParameter("Constant");
      par2[5] = gaus->GetParameter("Mean");
      par2[6] = gaus->GetParameter("Sigma");

      h_c_vector_[i]->GetXaxis()->SetRange(0, 0); // reset range

      // Temporarely subtract the pedestal gaussian to better fit the langaus

      vector<int> *reallocBin = new vector<int>;
      vector<double> *reallocBinContent = new vector<double>;

      for (int k = infbin; k < midbin + 20; k++)
      {
        double adjValue = h_c_vector_[i]->GetBinContent(k) - gaus->Eval(h_c_vector_[i]->GetBinCenter(k));
        reallocBin->push_back(k);
        reallocBinContent->push_back(h_c_vector_[i]->GetBinContent(k));
        if (adjValue >= 0)
        {
          h_c_vector_[i]->SetBinContent(k, adjValue);
        }
        else
        {
          h_c_vector_[i]->SetBinContent(k, 0.);
        }
      }

      // Evaluating the best starting parameters for the langaus

      par[2] = h_c_vector_[i]->Integral(); // Normalization factor

      j = h_c_vector_[i]->GetMaximumBin();
      par[1] = h_c_vector_[i]->GetBinCenter(j); // Most Probable value

      par[0] = h_c_vector_[i]->GetStdDev() * 0.1; // Width

      par[3] = 0.01; // Gaussian resolution

      // Fit of langaus function
      langaus->SetNpx(1000);
      langaus->SetParameters(par);
      langaus->SetParLimits(3, 0.001, 10);
      langaus->SetParNames("Width", "MP", "Area", "GSigma");

      h_c_vector_[i]->Fit(langaus, "Q", "same", -0.01, 0.05);

      // Saving the langaus fit parameters
      par2[0] = langaus->GetParameter("Width");
      par2[1] = langaus->GetParameter("MP");
      par2[2] = langaus->GetParameter("Area");
      par2[3] = langaus->GetParameter("GSigma");

      // Returning the histogram to normal
      for (int k = 0; k < reallocBin->size(); k++)
      {
        h_c_vector_[i]->SetBinContent((*reallocBin)[k], (*reallocBinContent)[k]);
      }
      delete reallocBin;
      delete reallocBinContent;

      // Fit of total distribution using the previous fit parameters (pedestal gaus + langaus)
      langausgaus->SetNpx(1000);
      langausgaus->SetParameters(par2);
      langausgaus->SetParLimits(3, 0.001, 10);
      langausgaus->SetParNames("Width", "MP", "Area", "G Sigma", "Ped Constant", "Ped Mean", "Ped Sigma");

      fitResults[i] = h_c_vector_[i]->Fit(langausgaus, "QS", "same", -0.01, 0.1);

      if ((int)fitResults[i] != 0)
      {
        fitError++;
        cout << "Fit of channel " << i << " was unsuccessful! Exit status " << (int)fitResults[i] << endl;
      }

      delete gaus;
      delete langaus;
      delete langausgaus;
    };
    delete c1;
    cout << "\r" << NCHANNELS - fitError << "/" << NCHANNELS << " fits were successful" << endl;
    return fitResults; // returns the number of unsuccessfull fits
  }
  else
  {
    cout << "Can't fit charge histograms." << endl;
    return nullptr;
  }
}

double AnaTools::lanGausFun(double *x, double *par)
{ // copied from laungaus.C example in Root documentation

  // Fit parameters:
  // par[0]=Width (scale) parameter of Landau density
  // par[1]=Most Probable (MP, location) parameter of Landau density
  // par[2]=Total area (integral -inf to inf, normalization constant)
  // par[3]=Width (sigma) of convoluted Gaussian function
  //
  // In the Landau distribution (represented by the CERNLIB approximation),
  // the maximum is located at x=-0.22278298 with the location parameter=0.
  // This shift is corrected within this function, so that the actual
  // maximum is identical to the MP parameter.

  // Numeric constants
  double invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
  double mpshift = -0.22278298;      // Landau maximum location

  // Control constants
  double np = 1000.0; // number of convolution steps
  double sc = 5.0;   // convolution extends to +-sc Gaussian sigmas

  // Variables
  double xx;
  double relsigma = 0;
  double mpc;
  double fland;
  double sum = 0.0;
  double xlow, xupp;
  double step;
  double i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Relative value of Gaussian sigma
  //relsigma = par[3]* (TMath::Sqrt(stepFun(x)) + 0.0001);

  // Range of convolution integral
  xlow = -0.01;
  xupp = 0.05;

  step = (xupp - xlow) / np;

  relsigma = par[3] * (TMath::Sqrt(TMath::Abs(x[0])));

  // Convolution integral of Landau and Gaussian by sum
  for (i = 1.0; i <= np / 2; i++)
  {
    xx = xlow + (i - .5) * step;
    
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(xx, x[0], relsigma);

    xx = xupp - (i - .5) * step;
    relsigma = par[3] * (TMath::Abs(x[0]));
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(xx, x[0], relsigma);
  }
  return (par[2] * step * sum * invsq2pi / relsigma);
}

double AnaTools::lanGausPlusGausFun(double *x, double *par)
{
  double invsq2pi = 0.3989422804014;
  double lanpar[4] = {par[0], par[1], par[2], par[3]};
  return lanGausFun(x, lanpar) + par[4] * TMath::Gaus(x[0], par[5], par[6])*invsq2pi/par[6];
}

double *AnaTools::EvaluateEfficiency()
{
  int goodFits = 0;
  static double result[NCHANNELS * 2] = {0};
  TFitResultPtr *fitResults = FitCharge();

  for (int i = 0; i < NCHANNELS; i++)
  {
    if ((int)fitResults[i] == 0)
      goodFits++;
  }

  if (goodFits != NCHANNELS)
  {
    cout << "Not all fits were successful! Can't evaluate channel efficiency." << endl;
    return nullptr;
  }

  for (int i = 0; i < NCHANNELS; i++)
  {
    double pedArea, langausArea, sigmaPed, sigmaLangaus;
    double *param;
    // double* covMatrix;
    // retrieve the the different fit functions
    TF1 *langausgaus = h_c_vector_[i]->GetFunction("langausgaus");
    TF1 *pedgaus = new TF1("pedgaus", "[0]*TMath::Gaus(x,[1],[2])", -2, 2);
    TF1 *langaus = new TF1("langaus", lanGausFun, -2, 2, 4);

    // Get fit parameters

    param = langausgaus->GetParameters();
    // covMatrix = fitResults[i]->GetCovarianceMatrix().GetMatrixArray();

    TMatrixDSym covMatrix = fitResults[i]->GetCovarianceMatrix();

    // Set the final fit parameters
    langaus->SetParameters(param);

    pedgaus->SetParameters(param + 4);

    // Calculate the area under the two functions
    pedArea = pedgaus->Integral(-0.01, 0.05);
    sigmaPed = pedgaus->IntegralError(-0.01, 0.05, param + 4, covMatrix.GetSub(4, 6, 4, 6).GetMatrixArray());

    langausArea = langaus->Integral(-0.01, 0.05);
    sigmaLangaus = langaus->IntegralError(-0.01, 0.05, param, covMatrix.GetSub(0, 3, 0, 3).GetMatrixArray());

    result[i] = langausArea / (langausArea + pedArea);
    result[i + NCHANNELS] = 1 / TMath::Power((langausArea + pedArea), 2) * TMath::Sqrt(TMath::Power(pedArea, 2) * TMath::Power(sigmaLangaus, 2) + TMath::Power(langausArea, 2) * TMath::Power(sigmaPed, 2));

    delete pedgaus;
    delete langaus;
  }

  return result;
}

void AnaTools::LoadInfo(TString infoFile)
{

  TFile *info = new TFile(infoFile, "READ");
  info->cd();

  TH1D *cutoff_values = (TH1D *)gDirectory->Get("cutoff_values");
  for (int i = 0; i < NCHANNELS; i++)
  {
    cutoff_[i] = cutoff_values->GetBinContent(i + 1);
  }

  info->Close();
  delete info;

  cout << "Loaded detector info from file!" << endl;

  outfile_->cd();

  return;
}

void AnaTools::SaveInfo(TString infoFile, TString infoType, TString mode){
  /*infoType options:
    - "cutoff";
    - "efficiency";
    - "gain".
  */

  if(infoFile == ""){
    infoFile = infoFile_;
  }

  // da fare controlla se il file già ci sta e in caso crea o modifica

  TFile *file = new TFile(infoFile, mode);
  file->cd();

  //Save cutoff values
  if(infoType.Contains("cutoff")){
    TH1D *h_cutoff = new TH1D("cutoff_values", "Cutoff Values", NCHANNELS, -0.5, NCHANNELS - 0.5);
    for (int i = 0; i < NCHANNELS; i++)
    {
      h_cutoff->SetBinContent(i + 1, cutoff_[i]);
    }

    file->Write();

    delete h_cutoff;
  }

  //Save efficiency values
  if(infoType.Contains("efficiency")){
    TH1D *h_efficiency = new TH1D("efficiency_values", "Efficiency Values", NCHANNELS, -0.5, NCHANNELS - 0.5);
    for (int i = 0; i < NCHANNELS; i++)
    {
      h_efficiency->SetBinContent(i + 1, efficiency_[i]);
    }

    file->Write();

    delete h_efficiency;
  }

  //Save gain values
  if(infoType.Contains("gain")){
    TH1D *h_gain = new TH1D("gain_values", "Gain Values", NCHANNELS, -0.5, NCHANNELS - 0.5);
    for (int i = 0; i < NCHANNELS; i++)
    {
      h_gain->SetBinContent(i + 1, gain_[i]);
    }

    file->Write();

    delete h_gain;
  }

  file->Close();
  delete file;

  return;
}

void AnaTools::SetCutoff(int i, const double &x)
{
  if (i < 0 || i >= NCHANNELS)
  {
    cout << "Can't set cutoff value, no " << i << "th channel exists." << endl;
    return;
  }

  cutoff_[i] = x;
  return;
}

void AnaTools::SetEfficiency(int i, const double &x){
  if (i < 0 || i >= NCHANNELS)
  {
    cout << "Can't set efficiency value, no " << i << "th channel exists." << endl;
    return;
  }

  efficiency_[i] = x;
  return;
}

void AnaTools::SetGain(int i, const double &x){
  if (i < 0 || i >= NCHANNELS)
  {
    cout << "Can't set egain value, no " << i << "th channel exists." << endl;
    return;
  }

  gain_[i] = x;
  return;
}

double *AnaTools::EvaluateMaxSignificanceBinCenter()
{
  if (bookings_[2] == 0)
  {
    cout << "Can't evaluate significance, charge histograms weren't booked." << endl;
    return nullptr;
  }

  // da fare: calcola fit se non è stato già fatto

  static double maxSignificanceCenter[NCHANNELS];

  for (int i = 0; i < NCHANNELS; i++)
  {

    TF1 *langausgaus = h_c_vector_[i]->GetFunction("langausgaus");
    TF1 *pedgaus = new TF1(Form("pedgaus_%i", i), "[0]*TMath::Gaus(x,[1],[2])", -0.01, 0.05);
    TF1 *langaus = new TF1(Form("langaus_%i", i), lanGausFun, -0.01, 0.05, 4);

    double *param = langausgaus->GetParameters();
    langaus->SetParameters(param);
    pedgaus->SetParameters(param + 4);

    TH1D *sig = new TH1D(Form("Significance_channel_%i", i), Form("Significance channel %i", i), 500, -0.01, 0.05);

    for (int j = 1; j <= 500; j++)
    {
      sig->SetBinContent(j, langaus->Integral(sig->GetBinCenter(j), 0.05, 1.0E-1) / TMath::Sqrt(langausgaus->Integral(sig->GetBinCenter(j), 0.05, 1.0E-1)));
    }

    maxSignificanceCenter[i] = sig->GetBinCenter(sig->GetMaximumBin());

    delete sig;
    delete pedgaus;
    delete langaus;
  }

  return maxSignificanceCenter;

}

void AnaTools::SetChargeRange(const double& xmin, const double& xmax){

  h_charge_range_[0] = xmin;
  h_charge_range_[1] = xmax;

  return;
}

double AnaTools::poisGausFun(double *x, double *par){
  /*
  There's a total of 6 parameters:
  - The pedestal expected value, mu_0;
  - The pedestal standard deviation, sigma_0;
  - The single photon charge expected value, mu_1;
  - The single photon charge standard deviation, sigma_1;
  - The poisson distribution expected value, lambda;
  - The integration constant.
  */

  double invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
 
  double coeff = par[5] * invsq2pi * TMath::Exp(-par[4]); 

  double sigma, sum = 0;

  for (int i = 0; i < 36; i++)
  {
    sigma = TMath::Sqrt(TMath::Power(par[1], 2) + TMath::Power(par[3], 2) * i);
    sum += TMath::Power(par[4], i) / (TMath::Factorial(i) * sigma) * TMath::Gaus(x[0], par[0]+i*par[2], sigma);
  }

  return sum * coeff;
}

double AnaTools::stepFun(double* x){ //step function
  if(x[0] >= 0.00001) return x[0];
  else return 0.00001;
} 

