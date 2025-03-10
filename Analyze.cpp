#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TSystem.h>
#include <stdio.h>
#include <TStyle.h>
#include <stdint.h>
#include "./Header/AnaTools.h"
#include "./Header/Event.h"

using namespace std;

int main(int argc, char *argv[]){

  string path("./");
  string inname("in.dat");
  string outname("out.root");

  int nevent=0;
  double cf=0.2;
  double th=0.;
  unsigned long int filepos=0;
  unsigned long int filelength;


  for(int i=0; i< argc; i++){
    if(strcmp("-path",argv[i])==0){ //Path for input and output file -- remember "/" at the end
      path.assign(argv[++i]);
    }
    if(strcmp("-in",argv[i])==0){ //input file name
      inname.assign(argv[++i]);
    }
    if(strcmp("-out",argv[i])==0){ //output file name
      outname.assign(argv[++i]);
    }
    if(strcmp("-nev",argv[i])==0){ //first event to analyze
      nevent = atoi(argv[++i]);
    }
     if(strcmp("-frac",argv[i])==0){ //value for constant fraction -- default 0.2
      cf = atof(argv[++i]);
    }
     if(strcmp("-th",argv[i])==0){ //value for threshold -- default 0.
      th = atof(argv[++i]);
    }
  }

  inname = path + inname;
  outname = path + outname;

  //open output file
  TFile *f = new TFile(outname.data(), "RECREATE"); //DISCLAIMER: RECREATE WILL OVERWRITE THE CONTENT OF THE OUTPUT FILE. IT'S BETTER TO DISTINGUISH EACH OUTPUT FILE BY THE CONSTANT OF TOF_cfm
  f->cd();

  //create an Event object
  Event * myEvent = new Event();


  //create an Analysis Tool object
  AnaTools *myAnaTools = new AnaTools(f,myEvent, cf, th, "info.root"); //info.root file contains some detector parameters

  //Calls for different kinds of analysis'
  myAnaTools->BookToF();
  myAnaTools->BookTime(); 
  myAnaTools->BookCharge();
  myAnaTools->BookWaveform();
  myAnaTools->BookPersistence();

  //open input file

  ifstream infile;
  infile.open(inname);
  if(!infile.is_open()){
    printf("*************** File %s not found, exiting from the program *****************\n",inname.data());
    exit(-1);
  }
  infile.seekg (0, infile.end);
  filelength = infile.tellg();

  
  //loop on events, read an event until file is finished
  while(filepos<filelength){
    myEvent->ReadEvent(inname, &filepos); 
    myAnaTools->Process(nevent);
    myEvent->Clear();
    nevent++;
    if(nevent%100==0){cout << "\r" <<"Processed " << 100*filepos/filelength << "%" << " of file" <<flush;}
  }

  cout << "\r" <<"Processed 100" << "%" << " of file" <<flush;
  cout << " " << endl;
	
  cout<< "Total No. Events read:" << nevent-1 << endl; 
 
  infile.close();
  
  //close output file
  f->Write();
  f->Close();

  delete myAnaTools;
  
  return 0;

}
