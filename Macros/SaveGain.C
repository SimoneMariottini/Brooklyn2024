#include "../Header/AnaTools.h"

void SaveGain(){
    AnaTools* MyAnaTools = new AnaTools("test_info.root"); //anatools object with path to info file
    
    //first set the gain values for every channel
    int channel = 1;
    double gain_value = 1;
    MyAnaTools->SetGain(0, );
    MyAnaTools->SetGain(1, );
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
    MyAnaTools->SetGain(15, );


    //then save info to file
    MyAnaTools->SaveInfo("gain");
}