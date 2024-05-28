#include "../Header/AnaTools.h"

void SaveGain(){
    AnaTools* MyAnaTools = new AnaTools("test_info.root"); //anatools object with path to info file
    
    //first set the gain values for every channel
    int channel = 1;
    double gain_value = 1;
    MyAnaTools->SetGain(channel, gain_value);

    //then save info to file
    MyAnaTools->SaveInfo("gain");
}