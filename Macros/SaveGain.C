#include "../Header/AnaTools.h"

#define NCHANNELS 16
#define elett 1.60217e-07 //carica dell'elettrone in picoCoulomb

void SaveGain(){
    AnaTools* MyAnaTools = new AnaTools("test_info.root"); //anatools object with path to info file
    
    //first set the gain values for every channel
    int channel = 1;
    double gain_value = 1;

    double charge[16];
    double gain[16];

    charge[0] = 1.794; //carica in picoCoulomb 
    charge[1] = 1.625;
    charge[2] = 1.57;
    charge[3] = 1.591;
    charge[4] = 1.603;
    charge[5] = 1.767;
    charge[6] = 1.802;
    charge[7] = 1.773;
    charge[8] = 1.756;
    charge[9] = 1.729;
    charge[10] = 1.691;
    charge[11] = 1.637;
    charge[12] = 1.635;
    charge[13] = 1.61;
    charge[14] = 1.628;
    charge[15] = 1.794;

    for (int i = 0; i < NCHANNELS; i++)
    {
        gain[i] = charge[i]/elett;
        MyAnaTools->SetGain(i, gain[i]);
    }
    
    //then save info to file
    MyAnaTools->SaveInfo("gain");
}




/* MyAnaTools->SetGain(1, );
    MyAnaTools->SetGain(2, );
    MyAnaTools->SetGain(3, );
    MyAnaTools->SetGain(4, );
    MyAnaTools->SetGain(5, );
    MyAnaTools->SetGain(6, );
    MyAnaTools->SetGain(7, );
    MyAnaTools->SetGain(8, );
    MyAnaTools->SetGain(9, );
    MyAnaTools->SetGain(10, );
    MyAnaTools->SetGain(11, );
    MyAnaTools->SetGain(12, );
    MyAnaTools->SetGain(13, );
    MyAnaTools->SetGain(14, );
    MyAnaTools->SetGain(15, );*/