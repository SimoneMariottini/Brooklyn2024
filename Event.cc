#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "Event.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <TMinuit.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <vector>
//#define NSAMPLING 1024

using namespace std;

Event::Event(){
 
  tot_charge = -1000;
  time = -1000;   
}

Event::~Event(){}


void Event::ReadEvent(string inname,unsigned long int* pos){

  ifstream myReadFile; //File "inname.dat" stream
  string fileline; //A line of the File
  string temp;    //Single value
  stringstream stream;    //Stream to navigate fileline
  vector<string> v;   //Vector of strings containing all the values (aka, temp)
  double va[NSAMPLING]; //Vector containing the values converted to double
  double vt[NSAMPLING]; //Vector containing the time steps
  int channel=0;
  //Reading position of the file
  myReadFile.open(inname);
  myReadFile.seekg(*pos);

  while(channel<16){  //Read file
    getline(myReadFile,  fileline); //write the line from file to string fileline
    if(fileline.find('=')==0){  //ignore lines which do not contain amplitude data
      fileline.clear();
      channel--;
    }
    else{
      stringstream stream(fileline);
      while (getline(stream, temp, ' ')) {    //Values in line are separated by space (' ')
	//	cout << fileline << endl;
	v.push_back(temp); //add value (still as string) to vector<string> temp
      }
      if(v.size()==NSAMPLING+1){
	v.pop_back();

      }//Sometimes, the program will a read an empty entry at the end. Delete that entry

      for (unsigned int i = 0; i < v.size(); i++) { //Convert vector<string> to array of double
        va[i]=stod(v[i]);
        vt[i]=i*SAMPLINGPERIOD;
      }
      Waveform *newchannel = new Waveform(&vt[0],&va[0]);
      myWaveforms.push_back(newchannel);
      v.clear();
    }

    channel++; 
  }

  *pos = myReadFile.tellg();
  //Clears and file closing, to prevent memory leaks
  v.clear();
  fileline.clear();
  myReadFile.close();
}


void Event::Clear(){


  for(long unsigned int iW=0;iW<myWaveforms.size();iW++){
    delete myWaveforms.at(iW);
  }
  myWaveforms.clear();

  tot_charge = -1000;
  time = -1000;

  return;

}



//GETTERS
vector<Waveform*> Event::getWaveforms(){
       return myWaveforms;
} 
double Event::GetTot_charge(){
       return tot_charge;
}
double Event::gettime(){
       return time;
}
 
//SETTERS
void Event::settot_charge(double caricaTotale){
       tot_charge=caricaTotale;
}
